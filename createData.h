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
using Tensor5 = vector<Tensor4>;

void Simple_Calculations(string incr, string folderName, double max_ratio, int minTradingDays, Matrix Rs)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "Simple";
    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor,
            iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);

    //Calculate files
    Matrix DataSet(0);
    for (int Period = 0; Period < iPeriods.size(); Period++)
    {
        //Calculating data
        Matrix beta_alpha_return = Beta_Alpha_Calculate(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods[Period], Inflation_Factor[Period], minTradingDays, false);
        push_back(DataSet, beta_alpha_return);
    }

    //Saving CSV and logfile
    SaveToCSV_transposed(folderName+"/"+methodName+"_CSV.txt", DataSet, fileNames);
    save_logfile(start, logMessage, folderName);

    //Saving individual files
    //for (int file = 0; file < fileNames.size(); ++file)   Save(folderName + fileNames[file], DataSet[file]);
}
void PrePost_Calculations(string incr, string folderName, double max_ratio, Matrix Rs)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "PrePost";
    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    int minTradingDays=0;
    Load_Data(Rs, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor,
            iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);

    //Calculate output data in twoPeriod_Data
    vector<Matrix> fullData(vector<Matrix>(2,Matrix(0,Vector(0))));
    for (int Period = 1; Period < iPeriods.size(); Period++)
    {
        //Calculating pre and post data
        vector<Matrix> prePostData = TwoPeriod_Calc(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods[Period-1], iPeriods[Period], {Inflation_Factor[Period-1], Inflation_Factor[Period]});
        push_back(fullData, prePostData);
    }

    //Saving CSV and logfile
    vector<string> prePost_str = {"_pre", "_post"};
    Save_TwoDataSet_CSV_transposed(folderName+"/"+methodName+"_CSV.txt", fullData, fileNames, prePost_str);
    save_logfile(start, logMessage, folderName);

    //Crating and saving individual dirs and files for backtesting use
    string preDirName = folderName + "/Pre_Period";
    string periodDirName = folderName + "/Period";
    vector<string> paths {preDirName, periodDirName};
    for (int prePost = 0; prePost < paths.size(); ++prePost) {
        mkdir(paths[prePost].c_str());
        for (int file = 0; file < fileNames.size(); ++file)
            Save(paths[prePost] + fileNames[file], fullData[prePost][file]);
        SaveToCSV_transposed(paths[prePost]+"/CSV" + prePost_str[prePost] + ".txt", fullData[prePost], fileNames);
    }
}
