

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
using Tensor5 [[maybe_unused]] = vector<Tensor4>;


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
        Data = Calculate_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[period], Inflation_Factor[period], 130, false);

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
            vector<Matrix> twoPeriod_Data = PrePost_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr], {Inflation_Factor[periodNr-1], Inflation_Factor[periodNr]});

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

 void old_BackTest(string incr, string folderName, string dataFolderName, Matrix Rs)
{
    //todo: Creating the BAB portfolio as if it was a stock (ETF) with daily returns
    auto start = std::chrono::system_clock::now();

    //Loading data
    string methodName = "BackTest_run";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);

    //Load data from prePost method
    vector <string> infoFiles = {"/year.txt", "/beta.txt", "/PERMNO.txt"};
    Tensor3 Data = Load_prePost(incr+"/"+dataFolderName, infoFiles);

    //Initiate readability variables
    int year = 0, beta = 1, PERMNO = 2;
    int pre = 0, post = 0;
    int low = 0, high = 1, all=2;

    //Low and high beta cutoff for portfolios
    double betaLow = 1, betaHigh = 1;

    //Initiate data variables
    int n_LHA = 3;
    int yrMin = Data[post][year][0];
    int yrMax = Data[post][year][Data[post][year].size()-1];
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(3, Intor(0)));
    Matrix beta_LHA(yrMax-yrMin+1, Vector(n_LHA,0));    //Low, High and All beta portfolio
    Vector (yrMax-yrMin+1);


    //Calculate which stocks should go in each (low/high beta) portfolio
    int st = 0;
    for (int yr = yrMin; yr <= yrMax; ++yr)
    {
        int iyr = yr-yrMin;
        while(Data[post][year][st] == yr)
        {
            //calculate beta
            double hist_Beta = Data[pre][beta][st];

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(Data[pre][PERMNO][st]);
            beta_LHA[iyr][all] += est_period_beta(hist_Beta);

            //Low beta portefolio
            if(hist_Beta < betaLow)
            {
                ID_LHA[iyr][low].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][low] += est_period_beta(hist_Beta);
            }
                //High beta portefolio
            else if(hist_Beta > betaHigh)
            {
                ID_LHA[iyr][high].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][high] += est_period_beta(hist_Beta);
            }
            st++;   //next stock
        }
        //Average beta in each portefolio
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            beta_LHA[iyr][LHA] /= ID_LHA[iyr][LHA].size();
    }

    Matrix LHA_return_beta_date(2*n_LHA+1, Vector(0));
    //Creating portfolios
    for (int iyr = 1; iyr <= yrMax-yrMin+1; ++iyr)
    {
        //Calculating period length
        int periodLength = iPeriods[iyr][1]-iPeriods[iyr][0]+1;

        for (int LHA = 0; LHA<2*n_LHA+1; ++LHA)
        {
            if(LHA < n_LHA)             //Calculating returns
                push_back(LHA_return_beta_date[LHA], Portfolio_Returns(ID_LHA[iyr-1][LHA], iDates, iPeriods[iyr], Rs));

            else if(2*n_LHA > LHA) {    //Calculating betas
                Vector v = Vector(periodLength, beta_LHA[iyr-1][LHA-n_LHA]);    //todo: make 1 line
                push_back(LHA_return_beta_date[LHA], v);
            }
            else                        //Calculating dates
                push_back(LHA_return_beta_date[LHA], int_to_double(dateRange(iPeriods[iyr], Dates)));
        }
    }

    //Saving CSV
    vector<string> headers = {"lowReturn", "highReturn", "allReturn", "lowPortBeta", "highPortBeta", "allPortBeta", "dates"};
    LHA_return_beta_date.push_back(skip_First_X(sp500, iPeriods[0][1]));
    LHA_return_beta_date.push_back(skip_First_X(riskFree, iPeriods[0][1]));
    headers.emplace_back("sp500");
    headers.emplace_back("riskFree");

    Save_CSV(folderName+"/BackTest_CSV.txt", LHA_return_beta_date, headers);
    save_logfile(start, logMessage, folderName);
}

 Vector Portfolio_Returns(Intor PERMNO, Intrix iDates, Intor iPeriod, Matrix Rs)
{
    //Initiating data structures
    int N = PERMNO.size();
    Vector iID(N);

    Intor iStart(N);
    Intor iRun_Period(N);

    int periodLength = iPeriod[1] - iPeriod[0] + 1;
    Vector portfolio_return(periodLength);
    Vector portfolio_kurs(periodLength);
    Vector portfolio_sum(periodLength);

    Vector shareprice(N,1);

    //Finding the indexes of the PERMNO's in Rs, to look stock up in Rs
    int ID_count = 0;
    for (int i = 0; i < Rs.size(); ++i)
        if(Rs[i][0] == PERMNO[ID_count])
        {
            iID[ID_count] = i;
            ID_count++;
        }

    //Finding the start and amount of trading days from start with stock data
    for (int i = 0; i < N; ++i)
    {
        int ID = iID[i];
        iStart  [i] = (iPeriod[0] - iDates[ID][1]) + 1;   //index where period start in Rs[i]
        iRun_Period[i] = Rs[ID].size() - iStart[i];   //index where period end in Rs[i]
    }

    //Calculating portfolio returns for period
    for (int day = 0; day<periodLength; ++day)
    {
        Vector shareprice_day(0);
        for (int stock = 0; stock<N; ++stock)
        {
            //If stock is still 'alive'
            if(day < iRun_Period[stock])
            {
                //Calculating new share price of portfolio
                double stock_return = Rs[iID[stock]][iStart[stock] + day];
                if(stock_return > -1.5) shareprice[stock] *= (1.0 + stock_return);
                shareprice_day.push_back(shareprice[stock]);
            }
        }
        //Daily portfolio shareprice and return is calculated

        //portfolio_kurs[day] = average(shareprice_day);  //why not shareprice, thus all stocks
        portfolio_kurs[day] = average(shareprice) * N/shareprice_day.size();

        if(day != 0)    portfolio_return[day] = portfolio_kurs[day]/portfolio_kurs[day-1] -1;
        else            portfolio_return[day] = portfolio_kurs[day]/1.0 - 1;
    }
    return portfolio_return;
}


 */

//old or unused calculate code
double CovSP500_log(Vector Stock, Vector sp500, int iStart_sp500, int iEnd_sp500, int iStart_DR)
{
    int d_iDates = iEnd_sp500 - iStart_sp500 + 1;
    double Mean_Stock = 0, Mean_sp500 = 0, Covv = 0;
    int trading_Days = 0;

    //Mean
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i+iStart_DR] != -2)
        {
            trading_Days++;
            Mean_Stock += log(1+Stock[i + iStart_DR]);
            Mean_sp500 += log(1+sp500[i + iStart_sp500]);
        }
    }
    Mean_Stock /= trading_Days;
    Mean_sp500 /= trading_Days;

    //CovSP500
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i+iStart_DR] != -2)
            Covv += (log(1+Stock[i+iStart_DR]) - Mean_Stock) * (log(1+sp500[i + iStart_sp500]) - Mean_sp500);
    }
    return Covv /(trading_Days - 1);
}
double Var_log(Vector sp500, int Start_number, int End_number)
{
    double Mean = 0, Variance = 0;
    //int d_dates = End_number-Start_number+1;
    //if(sp500.size() < d_dates) d_dates = sp500.size();

    int trading_Days=0;

    //Mean
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
        {
            trading_Days++;
            Mean += log(1+sp500[i]);
        }
    }
    Mean /= trading_Days;

    //Variance
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
            Variance += (log(1+sp500[i]) - Mean) * (log(1+sp500[i]) - Mean);
    }
    return Variance /(trading_Days - 1);
}
Intor Column_of_Intrix(Intrix A, int n)
{
    Intor v(A.size());
    for (int i = 0; i < A.size(); ++i) {
        v[i] = A[i][n];
    }
    return v;
}
bool areFilesExistInEitherDirectory(const vector<string>& filenames, const vector<string>& directoryPaths)
{
    bool allExist = true;
    for (const auto& filename : filenames) {
        for(const auto& directoryPath : directoryPaths)
        {
            const auto filepath = directoryPath + filename;
            if (!filesystem::exists(filepath)) {
                std::cerr << "File " << filepath << " does not exist\n";
                allExist = false;
            }
        }
    }
    if(allExist)    return true;
    else            return false;
}
Intrix dPeriods_From_iPeriod(Intrix iPeriod)
{
    Intor Dates = Load_Intor("Dly_DateList.txt");
    for(auto & i : iPeriod)
    {
        i[0] = Dates[i[0]];
        i[1] = Dates[i[1]];
    }
    return iPeriod;
}
int numFilesInDir(string dirPath)
{
    const std::filesystem::path dir_path = dirPath;
    int num_files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.path().filename().string()[0] != '.') {
            ++num_files;
        }
    }
    return num_files;
}
int countOfElements(Matrix A)
{
    int count = 0;
    for(auto v:A) {
        count += v.size();
    }
    return count;
}
Vector Portfolio_Returns_methods(Intor PERMNO, Intrix iDates, Intor iPeriod, Matrix Rs, Vector riskFree, string Method)
{
    //Initiating data structures
    int N = PERMNO.size();
    Vector iID(N);

    Intor iStart(N);
    Intor iRun_Period(N);

    int periodLength = iPeriod[1] - iPeriod[0] + 1;
    Vector portfolio_return(periodLength);
    Vector portfolio_sum(periodLength);

    Vector shareprice(N,1);

    //Finding the indexes of the PERMNO's in Rs, to look stock up in Rs
    int ID_count = 0;
    for (int i = 0; i < Rs.size(); ++i)
        if(Rs[i][0] == PERMNO[ID_count])
        {
            iID[ID_count] = i;
            ID_count++;
        }

    //Finding the start and amount of trading days from start with stock data
    for (int i = 0; i < N; ++i)
    {
        int ID = iID[i];
        iStart  [i] = (iPeriod[0] - iDates[ID][1]) + 1;   //index where period start in Rs[i]
        iRun_Period[i] = Rs[ID].size() - iStart[i];   //index where period end in Rs[i]
    }

    //Calculating portfolio returns for period
    for (int day = 0; day<periodLength; ++day)
    {
        //shareprice_last is updated
        Vector shareprice_last = shareprice;

        //shareprice is calculated
        for (int stock = 0; stock<N; ++stock)
        {
            //If stock is still 'alive', return is from stock, otherwise return is from the risk free rate
            if(day < iRun_Period[stock])    shareprice[stock] *= (1.0 + Rs[iID[stock]][iStart[stock] +day]);
            else if(Method == "RF") shareprice[stock] *= (1.0 + riskFree[iPeriod[0]+day]);
        }

        //Daily portfolio return is calculated based on either riskFree or EqualSplit when stock is missing data
        if(Method == "RF")
            portfolio_return[day] = sum(shareprice)/sum(shareprice_last) -1;
        else if(Method == "EQ")
        {
            double portPrice=0, portPrice_last=0;
            for (int stock = 0; stock<N; ++stock)
            {
                if(day < iRun_Period[stock])
                {
                    portPrice += shareprice[stock];
                    portPrice_last += shareprice_last[stock];
                }
            }
            portfolio_return[day] = portPrice/portPrice_last -1;
        }
        else if(Method == "Rebalance" || Method == "RB"){

            portfolio_return[day] = sum(shareprice)/sum(shareprice_last) -1;

            double sum_lost=0, sum_now=0, sum_last=0;
            for (int stock = 0; stock<N; ++stock)
            {
                if(day-1 < iRun_Period[stock])
                    sum_last += shareprice_last[stock];
                if(day < iRun_Period[stock])
                    sum_now += shareprice[stock];
                if(day-1 == iRun_Period[stock])
                    sum_lost += shareprice[stock];
            }
            double multiplier = sum_now/(sum_now-sum_lost);
            for (int stock = 0; stock<N; ++stock)
                if(day+1 < iRun_Period[stock])
                    shareprice[stock] *= multiplier;
        }
    }
    return portfolio_return;
}
void push_back(Intor& first, Intor last)
{
    first.insert(first.end(), last.begin(), last.end());
}
int day_after1926(int date)
{
    int year_after = date/10000 - 1926;
    int month_after = (date%10000)/100 - 1;
    int days_after = date%100-2;
    return 263 * year_after + month_after * 23 + days_after * 13/16;
} //TODO: Check if it works decently for years after 1927
/*int Slow_Find_iDate(int date, Intor Date_list)
{
    int index = 0;
    while(Date_list[index] < date)    index += 1;
    return index;
}*/


inline bool filePath_exists(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}
void Process_Files(string Dly_Mly_Both, bool Reload_Everything, bool Reload_smallFiles){
    bool Dly, Mly;  whatToLoad(Dly_Mly_Both, Dly, Mly);

    //Exo dir and required files
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(Dly_Mly_Both, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);


    vector<string> filenames = {"Dly_sp500.txt","Dly_DateList.txt", "Dly_Yly_RFR.txt", "Mly_pMC.txt", "First_MC.txt", "Yly_Inflation.txt"};
    if(Dly) filenames.push_back("Dly_Rs.txt");
    if(Mly) filenames.push_back("Mly_Rs.txt");

    if (areFilesExistInDirectory(filenames, Exo_FilePath)) {
        Intrix Rs_Dates;


        //Make dirs
        mkdir(Proccessed_FilePath.c_str());
        mkdir(Proccessed_Dly.c_str());
        mkdir(Proccessed_Mly.c_str());
        mkdir(Proccessed_Yly.c_str());

        cout << "\nCreating all necessary files and putting them in sub-directorys inside the directory " << Proccessed_FilePath << ":\n\n";

        //What to create files for:
        int max = 99999999;

        if(Dly)
        {
            //Create Daily big files:
            if (Reload_Everything || !areFilesExistInDirectory({"Rs.txt"}, Proccessed_Dly)) {
                Matrix Rs = Load_DR(Exo_FilePath+"Dly_Rs.txt", max);  //TODO: uncomment
                Compress_DR(Proccessed_Dly+"Rs.txt", Rs);   //TODO: uncomment
                cout << "Created " << Proccessed_Dly << "Rs.txt\n";
            }
            if (Reload_Everything || !areFilesExistInDirectory({"Rs_Dates.txt"}, Proccessed_Dly)) {
                Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
                Save(Proccessed_Dly+"Rs_Dates.txt", Rs_Dates);
                cout << "Created " << Proccessed_Dly << "Rs_Dates.txt\n";
            }
            else {
                Rs_Dates = Load_Intrix(Proccessed_Dly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }
            //Create Dly small files
            vector<string> smallFiles = {"Rs_iDates.txt", "iPeriods.txt", "riskFreeReturn.txt"};
            if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(smallFiles, Proccessed_Dly))
            {
                //DateList
                Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
                Dly_DateList.push_back(55555555);
                for(int i=0; i<20; i++) Dly_DateList.push_back(999999999);
                Save(Proccessed_Dly + "DateList.txt", Dly_DateList);

                //Rs_iDates
                Intrix Rs_iDates = Dly_Dates_to_iDates(Rs_Dates, Dly_DateList, 1);
                Save(Proccessed_Dly + "Rs_iDates.txt", Rs_iDates);

                //iPeriods
                Intrix iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);
                Save(Proccessed_Dly + "iPeriods.txt", iPeriods);

                //Risk free return
                Vector Dly_YearlyRFR = Load_Vector(Exo_FilePath + "Dly_Yly_RFR.txt");
                Vector Dly_RFR = DailyYearly_to_DailyDaily_Return(Dly_YearlyRFR, iPeriods);
                Save(Proccessed_Dly + "riskFreeReturn.txt", Dly_RFR);

                //sp500
                Vector sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
                Save(Proccessed_Dly+"sp500.txt",sp500);

                cout << "Created " << Proccessed_Dly << "Rs_iDates.txt\n";
                cout << "Created " << Proccessed_Dly << "iPeriods.txt\n";
                cout << "Created " << Proccessed_Dly << "riskFreeReturn.txt\n";
            }
        }
        if(Mly)
        {
            //Create Montly big files:
            if (Reload_Everything || !areFilesExistInDirectory({"Rs.txt"}, Proccessed_Mly)) {
                Matrix Rs = Load_DR(Exo_FilePath+"Mly_Rs.txt", max);  //TODO: uncomment
                Compress_DR(Proccessed_Mly+"Rs.txt", Rs);   //TODO: uncomment
                cout << "Created " << Proccessed_Mly << "Rs.txt\n";
            }
            if (Reload_Everything || !areFilesExistInDirectory({"Rs_Dates.txt"}, Proccessed_Mly)) {
                Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Mly_Rs.txt");
                Save(Proccessed_Mly+"Rs_Dates.txt", Rs_Dates);
                cout << "Created " << Proccessed_Mly << "Rs_Dates.txt\n";
            }
            else {
                Rs_Dates = Load_Intrix(Proccessed_Mly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }
            //Create Mly small files
            vector<string> smallFiles = {"Rs_iDates.txt", "iPeriods.txt", "riskFreeReturn.txt"};
            if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(smallFiles, Proccessed_Mly))
            {
                //DateList
                Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
                Intor Mly_DateList = Dly_to_Mly_DateList(Dly_DateList);
                Mly_DateList.push_back(55555555);
                for(int i=0; i<20; i++) Mly_DateList.push_back(999999999);
                Save(Proccessed_Mly + "DateList.txt", Mly_DateList);

                //Rs_iDates
                Intrix Rs_iDates = Mly_Dates_to_iDates(Rs_Dates, 1);  //TODO: Exact EOM date to iEOM function
                Save(Proccessed_Mly + "Rs_iDates.txt", Rs_iDates);

                //Mly_iPeriods
                Intrix Mly_iPeriods = x_iPeriods("Yly", "Mly", Mly_DateList);   //0-11, 12-23 , ...
                Save(Proccessed_Mly + "iPeriods.txt", Mly_iPeriods);

                //iPeriods needed to create 'Risk free return'
                Intrix Dly_M_iPeriods = x_iPeriods("Mly", "Dly", Dly_DateList);
                Intrix Dly_Y_iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);

                //Risk free return
                Vector Dly_YearlyRFR = Load_Vector(Exo_FilePath + "Dly_Yly_RFR.txt");
                Vector Dly_RFR = DailyYearly_to_DailyDaily_Return(Dly_YearlyRFR, Dly_Y_iPeriods);
                Vector Mly_RFR = period_accumulate_of_Dly(Dly_RFR, Dly_M_iPeriods);
                Save(Proccessed_Mly + "riskFreeReturn.txt", Mly_RFR);

                //sp500
                Vector sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
                Vector Mly_sp500 = period_accumulate_of_Dly(sp500, Dly_M_iPeriods);
                Save(Proccessed_Mly+"sp500.txt",Mly_sp500);
            }
            cout << "Created " << Proccessed_Mly << "Rs_iDates.txt\n";
            cout << "Created " << Proccessed_Mly << "iPeriods.txt\n";
            cout << "Created " << Proccessed_Mly << "riskFreeReturn.txt\n";
        }
        //Yly files, always make these if missing or small/big reload
        vector<string> Yly_smallFiles = {"pMC.txt", "Inflation_Factor.txt"};
        if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(Yly_smallFiles, Proccessed_Yly))
        {
            //Rs_Dates is needed for correcting Market Cap's otherwise 2 missing stocks
            if (Dly_Mly_Both != "Dly" || Dly_Mly_Both != "dly")
            {
                if (!filePath_exists(Proccessed_Dly+"Rs_Dates.txt"))
                {
                    Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
                    Save(Proccessed_Dly+"Rs_Dates.txt", Rs_Dates);
                    cout << "Created " << Proccessed_Dly << "Rs_Dates.txt\n";
                }
                else    Rs_Dates = Load_Intrix(Proccessed_Dly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }

            //Market Cap
            int factor = 1; //todo: needs to be 1 atm (Compress_MC)
            Matrix Mly_pMC = Load_Mly_MarketCap(Exo_FilePath+"Mly_pMC.txt", factor);
            Matrix First_MC = Load_Mly_MarketCap(Exo_FilePath + "First_MC.txt", factor);
            Matrix pMC_2comb = combine_First_with_Mly_MC(Mly_pMC, First_MC);
            Matrix pMC_3comb = Missing_MC_to_Zero(pMC_2comb, Rs_Dates);
            //Compress_MC(Proccessed_Mly + "pMC.txt", pMC_3comb, 1);
            Matrix Yly_pMC = MarketCap_Mly_to_Yly(pMC_3comb);
            Compress_MC(Proccessed_Yly + "pMC.txt", Yly_pMC, 1);


            //Inflation factors
            Vector Inflation = Load_Vector(Exo_FilePath + "Yly_Inflation.txt");
            Vector Inflation_factors = Inflation_Factors_from_yrly_inf(Inflation);
            Save(Proccessed_Yly +"Inflation_Factor.txt", Inflation_factors);
        }
    }
    else
        HowToGetStarted();
}


pair<vector<Intrix>, Matrix> StockSelector_old(vector<Matrix> Data, Vector betaCond)
{
    //Initiate readability variables
    int year = 0, beta = 1, PERMNO = 2;
    int pre = 0, post = 1;
    int low = 0, high = 1, all=2;

    //Initiate data variables
    int n_LHA = 3;
    int yrMin = Data[post][year][0];
    int yrMax = Data[post][year][Data[post][year].size()-1];
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Matrix beta_LHA(yrMax-yrMin+1, Vector(n_LHA,0));    //Low, High and All beta portfolio

    int st = 0;
    for (int yr = yrMin; yr <= yrMax; ++yr)
    {
        int iyr = yr-yrMin;
        while(Data[post][year][st] == yr)
        {
            //calculate beta
            double hist_Beta = Data[pre][beta][st];

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(Data[pre][PERMNO][st]);
            beta_LHA[iyr][all] += (hist_Beta);

            //Low beta portefolio
            if(hist_Beta > betaCond[0] && hist_Beta < betaCond[1])
            {
                ID_LHA[iyr][low].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][low] += (hist_Beta);
            }
                //High beta portefolio
            else if(hist_Beta > betaCond[2] && hist_Beta < betaCond[3])
            {
                ID_LHA[iyr][high].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][high] += (hist_Beta);
            }
            st++;   //next stock
        }
        //Average beta in each portefolio
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            beta_LHA[iyr][LHA] /= ID_LHA[iyr][LHA].size();
    }
    return {ID_LHA, beta_LHA};
}
Vector PortfolioReturns_method_old(Intor PERMNO, Intrix iDates, Intor iPeriod, Matrix Rs, Vector riskFree, const string& method)
{
    //Initiating data structures
    size_t N = PERMNO.size();
    Intor iID(N);

    vector<size_t> iStart(N);
    vector<size_t> iRun_Period(N);


    int periodLength = iPeriod[1] - iPeriod[0] + 1;
    Vector portfolio_return(periodLength);
    Vector portfolio_sum(periodLength);

    Vector shareprice(N,1);

    //Finding the indexes of the PERMNO's in Rs, to look stock up in Rs
    int ID_count = 0;
    for (int i = 0; i < Rs.size(); ++i)
        if(Rs[i][0] == PERMNO[ID_count])
        {
            iID[ID_count] = i;
            ID_count++;
        }

    //Finding the start and amount of trading days from start with stock data
    for (int i = 0; i < N; ++i)
    {
        int ID = iID[i];
        iStart[i] = (iPeriod[0] - iDates[ID][1]) + 1;   //index where period start in Rs[i]
        iRun_Period[i] = Rs[ID].size() - iStart[i];   //index where period end in Rs[i]
    }

    //Calculating portfolio returns for period
    for (int day = 0; day<periodLength; ++day)
    {
        if(method == "rebalance" || method=="RB")
        {
            double total_return=0;
            int totalLeft=0;
            for (int stock = 0; stock<N; ++stock) {
                if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)
                {
                    total_return += Rs[iID[stock]][iStart[stock] +day];
                    totalLeft++;
                }
            }
            portfolio_return[day] = total_return/totalLeft;
            continue;
        }
        //shareprice_last is updated
        Vector shareprice_last = shareprice;

        //shareprice is calculated
        for (int stock = 0; stock<N; ++stock)
        {
            //if nothing else, return is zero for the day
            double stock_return = 0;

            //If stock is still 'alive', and return is non empty (-2) return is from stock, otherwise return is from the risk free rate if method is selected
            if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)  stock_return = Rs[iID[stock]][iStart[stock] + day];
            else if(method == "RF" || method == "riskfree")                             stock_return = riskFree[iPeriod[0]+day];

            //Shareprice is now updated
            shareprice[stock] *= (1 + stock_return);
        }

        //Daily portfolio return is calculated based on either riskFree or distributed method when stock is missing data
        if(method == "RF" || method == "riskfree")
        {
            portfolio_return[day] = sum(shareprice)/sum(shareprice_last) -1;
        }
        else if(method == "DT" || method == "distributed")
        {
            double portPrice=0, portPrice_last=0;
            for (int stock = 0; stock<N; ++stock)
            {
                if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)
                {
                    portPrice += shareprice[stock];
                    portPrice_last += shareprice_last[stock];
                }
            }
            portfolio_return[day] = portPrice/portPrice_last -1;
        }
    }
    return portfolio_return;
}
Matrix LHA_PortfolioCreator_old(vector<Intrix> ID_LHA, Matrix beta_LHA, Matrix Rs, Vector riskFree, Intrix iPeriods, Intrix iDates, Intor Dates, int post_yrMin, int post_yrMax, string method)
{
    //Creating portfolios //add new function
    int n_LHA = ID_LHA[0].size();
    int iyr_start = post_yrMin-1926;
    Matrix LHA_PortData(2*n_LHA+1, Vector(0));
    for (int iyr = iyr_start; iyr <= post_yrMax-1926; ++iyr)  //use yrTotal instead
    {
        //Calculating period length
        int periodLength = iPeriods[iyr][1]-iPeriods[iyr][0]+1;

        //Calculating returns
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA], PortfolioReturns_method_old(ID_LHA[iyr-iyr_start][LHA], iDates, iPeriods[iyr], Rs, riskFree, method));

        //Calculating betas
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA+n_LHA], Vector(periodLength, beta_LHA[iyr-iyr_start][LHA]));

        //Calculating dates
        push_back(LHA_PortData[2*n_LHA], int_to_double(dateRange(iPeriods[iyr], Dates)));
    }
    return LHA_PortData;
}


/*Intrix Dly_Dates_to_iDates(int (*Find_iDate)(int, Intor) ,Intrix Rs_Dates, Intor DateList, int start_j)
{
    for (auto & Date : Rs_Dates)
    {
        for (int j = start_j; j < Date.size(); ++j)
        {
            Date[j] = Find_iDate(Date[j], DateList);
        }
    }
    return Rs_Dates;
}*/
/*vector<Intrix> SplitPeriods(const Intrix& iPeriods, int x, bool rest_in_first_period)
{
    int numRowsPerVec = floor(iPeriods.size() / x);
    int rest = iPeriods.size() - x * numRowsPerVec;

    vector<vector<vector<int>>> threeDimVec(x, vector<vector<int>>(numRowsPerVec, vector<int>(2)));

    int currRow = 0;

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < numRowsPerVec; j++) {
            threeDimVec[i][j][0] = iPeriods[currRow][0];
            threeDimVec[i][j][1] = iPeriods[currRow][1];
            currRow++;
        }
        if(i==0)
        {
            if(rest_in_first_period)
            {
                if (rest > 0)
                {
                    for (int s = 0; s < rest; ++s)
                    {
                        threeDimVec[0].push_back({iPeriods[currRow][0],iPeriods[currRow][1]});
                        currRow++;
                    }
                }
            }
        }
    }

    cout << "Periods of years was made to be (" << x << " x " << numRowsPerVec << ")";
    if(rest > 0)
    {
        if(rest_in_first_period)
            cout << ", first period have " << numRowsPerVec+rest << " periods.\n";
        else
            cout << ", first period also have " << numRowsPerVec << " periods, so " << rest << " yeas was not included in any period.\n";
    }
    else
        cout << ".\n";

    return threeDimVec;
}*/


/*
Matrix LHA_PortData;

if(beta_weighted)
{
    //Beta weighted
    auto [ID_LHA, beta_LHA] = StockSelector_stockBeta(Data, betaCond);
    LHA_PortData = LHA_PortfolioCreator_stockBeta(ID_LHA, beta_LHA, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin, post_yrMax, method, beta_weighted);
}
else
{
    //Unweighted
    auto [ID_LHA, beta_LHA] = StockSelector(Data, betaCond);
    LHA_PortData = LHA_PortfolioCreator(ID_LHA, beta_LHA, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin, post_yrMax, method);
}
*/

//old method (unweighted)   //todo: remove
//auto [ID_LHA, beta_LHA] = StockSelector(Data, betaCond);
//Matrix LHA_PortData = LHA_PortfolioCreator(ID_LHA, beta_LHA, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin, post_yrMax, method);
