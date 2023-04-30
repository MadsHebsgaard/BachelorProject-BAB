

#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()


#include <chrono>
#include <ctime>


#pragma once

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;


/*

//createData
void Era_Calculations_test(string inc, string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    string methodName = "Era";

    vector<string> logMessage, fileNames;
    string Exo_FilePath, Proccessed_FilePath;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;

    Load_Data(DR, max_ratio, minTradingDays, logMessage, Exo_FilePath, Proccessed_FilePath, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, inc);
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
            Matrix beta_alpha_return = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr], Inflation_Factor[periodNr], minTradingDays, false);
            push_back(DataSet, beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(Era + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], DataSet[file]);
    }
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
}
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
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
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
            Matrix beta_alpha_return = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr], Inflation_Factor[periodNr], minTradingDays, false);
            push_back(DataSet, beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(Era + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], DataSet[file]);
    }
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
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
        Data = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[period], Inflation_Factor[period], 130, false);

        // Creating a directory and files for the period
        string periodDirName = dirName + "_" + to_string(period + 1 - skipFirst_n);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], Data[file]);
    }
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
}
void Era_Period_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR) {
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
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
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
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
}


//backtest
void BackTesting(string folderName, string dataFolderName)
{
    int iBeta = 0;
    int iReturn = 1;
    int iAlpha = 2;
    int iSP500 = 3;
    int iRiskFree = 4;

    vector <string> infoFiles = {"/beta.txt", "/akk_return.txt", "/alpha.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/PERMNO.txt"};
    Tensor5 Data = Load_Era_Period_PrePost(dataFolderName, infoFiles);

    double betaLow = 0.5;
    double betaHigh = 1.25;

    //todo: Max_high_Beta krav?
    //todo: Min_low_Beta krav?

    int stockInf = infoFiles.size() - 1; //Permno is excluded
    int iCount = 5;
    int eraCount = Data.size();
    int maxPeriodCount = Data[0].size();

    vector<Matrix> lowPortfolio(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));
    vector<Matrix> highPortfolio(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));
    vector<Matrix> BAB(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));

    vector<Matrix> preBeta(Data.size(), Matrix (maxPeriodCount, Vector(2, 0)));
    vector<Matrix> Count(Data.size(), Matrix (maxPeriodCount, Vector(2, 0)));

    Matrix BAB_Average(eraCount, Vector(stockInf+1, 0));



    //TODO: Create an Intor of PERMNOS to go in High and Low beta portfolios
    //TODO: Creating the BAB portfolio as if it was a stock (ETF) with daily returns
    //TODO: plot the BAB performance VS sp500 performance


    for (int eraNr = 0; eraNr < eraCount; ++eraNr) {
        int periodCount = Data[eraNr].size();   //Todo: removed -1
        Data[eraNr].resize(periodCount);

        for (int periodNr = 0; periodNr < periodCount; ++periodNr) {
            for (int stockNr = 0; stockNr < Data[eraNr][periodNr][0][0].size(); ++stockNr) {

                double hist_Beta = Data[eraNr][periodNr][0][0][stockNr];
                //Low beta portefolio
                if(hist_Beta < betaLow)
                {
                    double ammount = betaLow - hist_Beta;
                    ammount = 1;    //ammount = 1 for equal weighted
                    preBeta[eraNr][periodNr][0] += ammount * hist_Beta;
                    Count[eraNr][periodNr][0] += ammount;
                    for (int inform = 0; inform < stockInf; ++inform) {
                        lowPortfolio[eraNr][periodNr][inform] += ammount * Data[eraNr][periodNr][1][inform][stockNr];
                    }
                }
                    //High beta portefolio
                else if(hist_Beta > betaHigh)
                {
                    double ammount = hist_Beta - betaHigh;  //ammount = 1 for equal weighted
                    ammount = 1;    //ammount = 1 for equal weighted
                    preBeta[eraNr][periodNr][1] += ammount * hist_Beta;
                    Count[eraNr][periodNr][1] += ammount;
                    for (int inform = 0; inform < stockInf; ++inform) {
                        highPortfolio[eraNr][periodNr][inform] += ammount * Data[eraNr][periodNr][1][inform][stockNr];
                    }
                }
            }
            preBeta[eraNr][periodNr][0] /= Count[eraNr][periodNr][0];
            preBeta[eraNr][periodNr][1] /= Count[eraNr][periodNr][1];

            for (int info = 0; info < stockInf; ++info) {
                lowPortfolio[eraNr][periodNr][info] /= Count[eraNr][periodNr][0];     //Average
                highPortfolio[eraNr][periodNr][info] /= Count[eraNr][periodNr][1];    //Average
            }
            //double LowAmmount = preBeta[periodNr][1] / (preBeta[periodNr][0]+preBeta[periodNr][1]);
            double BhighLow = preBeta[eraNr][periodNr][1] - preBeta[eraNr][periodNr][0];
            double lowAmmount = ( preBeta[eraNr][periodNr][1] / BhighLow);
            double highAmmount = -( preBeta[eraNr][periodNr][0] / BhighLow);

            for (int info = 0; info < stockInf; ++info) {
                BAB[eraNr][periodNr][info] = lowAmmount * lowPortfolio[eraNr][periodNr][info] + highAmmount * highPortfolio[eraNr][periodNr][info];
                BAB_Average[eraNr][info] += BAB[eraNr][periodNr][info];
            }
            BAB_Average[eraNr][iCount] += Count[eraNr][periodNr][0] + Count[eraNr][periodNr][1];
        }
        for (int info = 0; info < stockInf+1; ++info)   BAB_Average[eraNr][info] /= periodCount;
    }



    //Save BAB, low and high portfolio
    mkdir("Data/BackTest");
    string folderPath = "Data/BackTest/" + folderName;
    std::filesystem::remove_all(folderPath.c_str());
    mkdir(folderPath.c_str());
    vector<string> backTestDirs {"/BAB portfolio", "/LowBeta portfolio", "/HighBeta portfolio"};

    for (int eraNr = 0; eraNr < eraCount; ++eraNr) {
        string eraDir = folderPath + "/Era_" + to_string(eraNr+1);
        mkdir(eraDir.c_str());

        vector<string> btDirs_copy = backTestDirs;
        for(auto& dir:btDirs_copy)
        {
            dir = eraDir + dir;
            mkdir(dir.c_str());
        }
        //BAB portfolio
        for (int info = 0; info < stockInf; ++info) {
            Matrix BAB_Data(3,Vector(Data[eraNr].size()));

            for (int periodNr = 0; periodNr < Data[eraNr].size(); ++periodNr) {
                BAB_Data[0][periodNr] =           BAB[eraNr][periodNr][info];
                BAB_Data[1][periodNr] =  lowPortfolio[eraNr][periodNr][info];
                BAB_Data[2][periodNr] = highPortfolio[eraNr][periodNr][info];
            }
            for (int dirNr = 0; dirNr < btDirs_copy.size(); ++dirNr) {
                string infoDir = btDirs_copy[dirNr] + infoFiles[info];
                Save(infoDir, BAB_Data[dirNr]);
            }
        }
    }
    //TODO: make log file
    return;
}


//Other
void SaveToCSV_DataSet(string fn, Matrix A, vector<string> headers)
{
    char delim = ';';
    ofstream fil(fn);
    int n = headers.size();
    for (int i = 0; i<n-1; ++i) {
        fil << headers[i] << delim;
    }
    fil << headers[n-1] << endl;

    for (int i = 0; i<A[0].size(); ++i) {
        for (int j = 0; j<n-1; ++j)
            fil << A[i][j] << delim;
        fil << A[i][n-1] << endl;
    }
}
void Load_Data(Matrix& DR, double& max_ratio, int& minTradingDays, vector<string>& logMessage
        , string& Exo_FilePath, string& Proccessed_FilePath, Matrix& MC, Vector& Inflation_Factor
        , Intrix& iDates, Vector& sp500, Vector& riskFree, Intor& Dates, Intrix& iPeriods
        , string& folderName, string& methodName, vector<string>& fileNames)
{
    //File paths
    Exo_FilePath = "Data/Input/Exo_Files/";
    Proccessed_FilePath = "Data/Input/Processed_Files/";
    fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor.txt", "/year.txt"};

    //Log messages
    logMessage = {};
    logMessage.push_back("max_ratio = "+to_string(max_ratio));
    logMessage.push_back("minTradingDays = "+formatNumber(minTradingDays));
    logMessage.push_back("DR.size() = "+formatNumber(DR.size()));
    string elementsDR = formatNumber(countOfElements(DR)-DR.size());

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);
    string elementsDR_ny = formatNumber(countOfElements(DR)-DR.size());
    logMessage.push_back("DR_ny.size() = "+formatNumber(DR.size()));


    //Market Cap
    MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR
    iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
    riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Dates = Load_Intor(Exo_FilePath+"Dly_DateList.txt");

    //Load iPeriods and create Era_List
    iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    logMessage.push_back("Data period: " + formatDate(Dates[iPeriods[0][0]]) + " - " + formatDate(Dates[iPeriods[iPeriods.size()-1][1]]));
    logMessage.push_back("Returns in DR = " + elementsDR);
    logMessage.push_back("Returns in DR_ny = " + elementsDR_ny);


    //Create files and folderName = folderPath
    mkdir("Data/Output");
    cout << methodName << "_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir(("Data/Output/" + methodName).c_str());
    folderName = "Data/Output/" + methodName + "/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());
}


 */