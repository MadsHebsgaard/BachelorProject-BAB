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

Matrix Simple_Calculations(Matrix Rs, Matrix MC, Intrix iDates, Vector sp500, Vector riskFree, Intor Dates, Intrix iPeriods, Vector Inflation_Factor, int minTradingDays)
{
    //Calculate DataSet
    Matrix DataSet(0);
    for (int Period = 0; Period < iPeriods.size(); Period++)
    {
        //Calculating data
        Matrix Yly_Data = Calculate_Performance(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods[Period],
                                                   Inflation_Factor[Period], minTradingDays, false);
        //Inserting Data
        push_back(DataSet, Yly_Data);
    }
    return DataSet;
}
void Simple_run(string incr, string folderName, int minTradingDays, Matrix Rs)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "Simple";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);

    //Calculate DataSet
    Matrix DataSet = Simple_Calculations(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods, Inflation_Factor, minTradingDays);

    //Saving CSV and logfile
    Save_CSV(folderName+"/"+methodName+"_CSV.txt", DataSet, fileNames);
    save_logfile(start, logMessage, folderName);

    //Saving individual files
    //for (int file = 0; file < fileNames.size(); ++file)   Save(folderName + fileNames[file], DataSet[file]);
}



vector<Matrix> PrePost_Calculations(Matrix Rs, Matrix MC, Intrix iDates, Vector sp500, Vector riskFree, Intor Dates, Intrix iPeriods, Vector Inflation_Factor)
{
    //Calculate twoPeriod_Data
    vector<Matrix> fullData(vector<Matrix>(2,Matrix(0,Vector(0))));
    for (int Period = 1; Period < iPeriods.size(); Period++)
    {
        //Calculating pre and post data
        vector<Matrix> Yly_Data = PrePost_Performance(Rs, MC, iDates, sp500, riskFree, Dates,
                {iPeriods[Period-1], iPeriods[Period]}, {Inflation_Factor[Period-1], Inflation_Factor[Period]});
        push_back(fullData, Yly_Data);
    }
    return fullData;
}
void PrePost_run(string incr, string folderName, Matrix Rs)
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
    vector<Matrix> fullData = PrePost_Calculations(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods, Inflation_Factor);


    //Saving CSV and logfile
    vector<string> prePost_str = {"_pre", "_post"};
    Save_TwoDataSet_CSV(folderName+"/"+methodName+"_CSV.txt", fullData, fileNames, prePost_str);
    save_logfile(start, logMessage, folderName);

    //Crating and saving individual dirs and files for backtesting use
    vector<string> paths {folderName + "/Pre_Period", folderName + "/Period"};
    for (int prePost = 0; prePost < paths.size(); ++prePost) {
        mkdir(paths[prePost].c_str());
        for (int file = 0; file < fileNames.size(); ++file)
            Save(paths[prePost] + fileNames[file], fullData[prePost][file]);
        Save_CSV(paths[prePost]+"/CSV"+prePost_str[prePost]+".txt", fullData[prePost], fileNames);
    }
}