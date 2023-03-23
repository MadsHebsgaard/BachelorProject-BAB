#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()


using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;

#pragma once

void Era_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    string  methodName = "Era_PrePost";


    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates
            , sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames);


    int minEraLength = floor((iPeriods.size()-1.0)/(n_Eras+0.0));
    int rest = iPeriods.size()-1-minEraLength*n_Eras;
    int eraLength;
    int periodNr = 1;


    for (int Era = 0; Era < n_Eras; ++Era)
    {
        vector<Matrix> tP_DataSet(vector<Matrix>(2,Matrix(0,Vector(0))));

        if(Era < rest)  eraLength = minEraLength+1;
        else            eraLength = minEraLength;

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 1; Period <= eraLength; Period++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            vector<Matrix> twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr], {Inflation_Factor[periodNr-1],Inflation_Factor[periodNr]});
            push_back(tP_DataSet, twoPeriod_Data);
            periodNr++;
        }
        //Create directory and path for era
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());
        string preDirName = dirName + "/Pre_Period";
        string periodDirName = dirName + "/Period";

        //Create files in directory for the era
        vector<string> paths {preDirName, periodDirName};
        for (int prePost = 0; prePost < paths.size(); ++prePost) {
            mkdir(paths[prePost].c_str());
            for (int file = 0; file < fileNames.size(); ++file)
                Save(paths[prePost] + fileNames[file], tP_DataSet[prePost][file]);
        }
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
void Era_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    string methodName = "Era";


    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames);
    string dirName = folderName + "/Era";


    int minEraLength = floor((iPeriods.size()-1.0)/(n_Eras+0.0));
    int rest = iPeriods.size()-1-minEraLength*n_Eras;
    int eraLength;
    int periodNr = 0;

    for (int Era = 0; Era < n_Eras; ++Era)
    {
        Matrix DataSet(0);
        if(Era < rest)  eraLength = minEraLength+1;
        else            eraLength = minEraLength;

        for (int Period = 0; Period < eraLength; Period++)
        {
            periodNr++;
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            Matrix beta_alpha_return = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr], Inflation_Factor[periodNr], minTradingDays);
            push_back(DataSet, beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(Era + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], DataSet[file]);
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
void Period_Calculations(string folderName, double max_ratio, int minTradingDays, Matrix DR, int skipFirst_n)
{
    string methodName = "Period";


    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates
            , sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames);

    string dirName = folderName + "/Period";
    Matrix Data;

    for (int period = skipFirst_n; period < iPeriods.size(); ++period)
    {
        //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
        Data = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[period], Inflation_Factor[period], 130);

        // Creating a directory and files for the period
        string periodDirName = dirName + "_" + to_string(period + 1 - skipFirst_n);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], Data[file]);
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}

void Era_Period_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)   //TODO: Make this
{
    string methodName = "Era_Period_PrePost";


    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates
            , sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames);

    int minEraLength = floor((iPeriods.size()-1.0)/(n_Eras+0.0));
    int rest = iPeriods.size()-1-minEraLength*n_Eras;
    int eraLength;
    int periodNr = 1;

    for (int Era = 0; Era < n_Eras; ++Era)
    {
        //Create directory and path for era
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());

        if(rest <= n_Eras){
            if(Era < rest)  eraLength = minEraLength+1;
            else            eraLength = minEraLength;
        }
        else if(Era == 0)   eraLength = minEraLength+rest;

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int eraPeriod = 1; eraPeriod <= eraLength; eraPeriod++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            vector<Matrix> twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr - 1], iPeriods[periodNr], {Inflation_Factor[periodNr-1], Inflation_Factor[periodNr]});

            //Create dirs and then files
            string preDirName = dirName + "/Period_" + to_string(eraPeriod);
            string periodDirName = dirName + "/Period_" + to_string(eraPeriod);
            mkdir(preDirName.c_str());
            mkdir(periodDirName.c_str());
            preDirName = preDirName + "/Pre_Period";
            periodDirName = periodDirName + "/Period";

            vector<string> prePost {preDirName, periodDirName};
            for (int i = 0; i < prePost.size(); ++i) {
                mkdir(prePost[i].c_str());
                for (int file = 0; file < fileNames.size(); ++file)
                    Save(prePost[i] + fileNames[file], twoPeriod_Data[i][file]);
            }
            periodNr++;
        }
        cout << "Era " << Era << " is done.\n";
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
void Era_PrePost_Period_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    string methodName = "Era_PrePost_Period";


    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates
            , sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames);

    int minEraLength = floor((iPeriods.size()-1.0)/(n_Eras+0.0));
    int rest = iPeriods.size()-1-minEraLength*n_Eras;
    int eraLength;
    int periodNr = 1;

    for (int Era = 0; Era < n_Eras; ++Era)
    {
        //Create directory and path for era + pre(period)
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());
        string preDirName = dirName + "/Pre_Period";
        string postDirName = dirName + "/Period";
        mkdir(preDirName.c_str());
        mkdir(postDirName.c_str());

        if(Era < rest)  eraLength = minEraLength+1;
        else            eraLength = minEraLength;

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int eraPeriod = 1; eraPeriod <= eraLength; eraPeriod++)
        {
            //if(Era == 0 && eraPeriod == 0) continue;
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            vector<Matrix> twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr], {Inflation_Factor[periodNr-1], Inflation_Factor[periodNr]});

            //Create dirs and then files
            string periodPre_DirName = preDirName + "/Period_" + to_string(eraPeriod);
            string periodPost_DirName = postDirName + "/Period_" + to_string(eraPeriod);

            mkdir(periodPre_DirName.c_str());
            mkdir(periodPost_DirName.c_str());

            vector<string> prePost {periodPre_DirName, periodPost_DirName};
            for (int i = 0; i < prePost.size(); ++i) {
                mkdir(prePost[i].c_str());
                for (int file = 0; file < fileNames.size(); ++file)
                    Save(prePost[i] + fileNames[file], twoPeriod_Data[i][file]);
            }
            periodNr++;
        }
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}