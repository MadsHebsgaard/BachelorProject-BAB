#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()
#pragma once

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;

Matrix Simple_Calculations(const Matrix& Rs, const Matrix& MC, const Intrix& iDates,
        const Vector& sp500, const Vector& riskFree, const Intor& Dates,
        Intrix iPeriods, Vector Inflation_Factor, int minTradingDays, bool logarithm)
{
    //Calculate DataSet
    Matrix DataSet(0);
    for (int Period = 0; Period < iPeriods.size(); Period++)
    {
        //Calculating data
        Matrix Yly_Data = Calculate_Performance(Rs, MC, iDates, sp500, riskFree, Dates,
                iPeriods[Period], Inflation_Factor[Period], minTradingDays, false, logarithm);
        //Inserting Data
        push_back(DataSet, Yly_Data);
    }
    return DataSet;
}
void Simple_run(const string& incr, string folderName, int minTradingDays, Matrix Rs, bool logarithm)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "Simple";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree,
            Dates, iPeriods, folderName, methodName, fileNames, incr);

    //Calculate DataSet
    Matrix DataSet = Simple_Calculations(Rs, MC, iDates, sp500, riskFree,
            Dates, iPeriods, Inflation_Factor, minTradingDays, logarithm);

    //Saving CSV and logfile
    Save_CSV(folderName+"/"+methodName+"_CSV.txt", DataSet, fileNames);
    save_logfile(start, logMessage, folderName, logarithm);
}


vector<Matrix> PrePost_Performance(const Matrix& DR, const Matrix& MC, const Intrix& iDates, const Vector& sp500,
        const Vector& riskFree, const Intor& Dates, Intrix iPeriod, const Vector& Inflation_factor, bool logarithm)
{
    vector<Matrix> A(0);
    //size of BetaOverlap is a tradeoff of bias vs robust, makes Post_beta more robust at the cost of making Post_beta biased towards Pre_beta. Tradeoff is probably not bad.
    int BetaPrePeriod = 500, BetaOverlap = 250;  //Dly
    if (iPeriod[1][0] - iPeriod[0][0] == 12)    //Mly
    {
        BetaPrePeriod = 24;
        BetaOverlap = 12;
    }

    Matrix Pre_Data = Calculate_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[0],
            Inflation_factor[0], BetaPrePeriod, true, logarithm);
    Matrix Post_Data = Calculate_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[1], Inflation_factor[1],
            BetaOverlap, false, logarithm);
    return Overlapping_ID_Matrix(Pre_Data, Post_Data, 3);
}
vector<Matrix> PrePost_Calculations(const Matrix& Rs, const Matrix& MC, const Intrix& iDates, const Vector& sp500,
        const Vector& riskFree, const Intor& Dates, const Intrix& iPeriods, const Vector& Inflation_Factor, bool logarithm)
{
    //Calculate twoPeriod_Data
    vector<Matrix> fullData(vector<Matrix>(2,Matrix(0,Vector(0))));
    for (int Period = 1; Period < iPeriods.size(); Period++)
    {
        //Calculating pre and post data
        vector<Matrix> Yly_Data = PrePost_Performance(Rs, MC, iDates, sp500, riskFree, Dates,
                {iPeriods[Period-1], iPeriods[Period]}, {Inflation_Factor[Period-1], Inflation_Factor[Period]}, logarithm);
        push_back(fullData, Yly_Data);
    }
    return fullData;
}
void PrePost_run(const string& incr, string folderName, Matrix Rs, bool logarithm)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "PrePost";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);


    //Calculate output data in twoPeriod_Data
    vector<Matrix> fullData = PrePost_Calculations(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods, Inflation_Factor, logarithm);


    //Saving CSV and logfile
    vector<string> prePost_str = {"_pre", "_post"};
    Save_TwoDataSet_CSV(folderName+"/"+methodName+"_CSV.txt", fullData, fileNames, prePost_str);
    save_logfile(start, logMessage, folderName, logarithm);

    //Crating and saving individual dirs and files for backtesting use
    vector<string> paths {folderName + "/Pre_Period", folderName + "/Period"};
    for (int prePost = 0; prePost < paths.size(); ++prePost) {
        mkdir(paths[prePost].c_str());
        for (int file = 0; file < fileNames.size(); ++file)
            Save(paths[prePost] + fileNames[file], fullData[prePost][file]);
        Save_CSV(paths[prePost]+"/CSV"+prePost_str[prePost]+".txt", fullData[prePost], fileNames);
    }
}
