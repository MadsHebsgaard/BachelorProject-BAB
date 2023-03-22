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
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);

    //Market Cap
    Matrix MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Vector Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);


    cout << "Era_PrePost_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Era_PrePost");
    folderName = "Data/Output/Era_PrePost/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;


    int minEraLength = floor((iPeriods.size()-1.0)/(n_Eras+0.0));
    int rest = iPeriods.size()-1-minEraLength*n_Eras;
    int eraLength;
    int periodNr = 1;

    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor"};

    for (int Era = 0; Era < n_Eras; ++Era)
    {
        vector<Matrix> tP_DataSet(vector<Matrix>(2,Matrix(0,Vector(0))));

        if(Era < rest)  eraLength = minEraLength+1;
        else            eraLength = minEraLength;

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 1; Period <= eraLength; Period++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr], {Inflation_Factor[periodNr-1],Inflation_Factor[periodNr]});
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
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);

    //Market Cap
    Matrix MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Vector Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Periods
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    //vector<Intrix> Era_List = SplitPeriods(iPeriods, n_Eras, true);

    cout << "Files was Loaded.\n\n";

    Matrix beta(n_Eras, Vector(0)), alpha(n_Eras, Vector(0)), PERMNO(n_Eras, Vector(0));
    Matrix akk_return(n_Eras, Vector(0)), akk_sp500(n_Eras, Vector(0)), akk_riskFree(n_Eras, Vector(0));
    Matrix beta_alpha_return;

    mkdir("Data/Output");
    mkdir("Data/Output/Era");
    folderName = "Data/Output/Era/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());

    string dirName = folderName + "/Era";
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor"};

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
            beta_alpha_return = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr], Inflation_Factor[periodNr], minTradingDays);
            push_back(DataSet, beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(Era + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file) Save(periodDirName + fileNames[file], DataSet[file]);
    }
}
void Period_Calculations(string folderName, double max_ratio, int minTradingDays, Matrix DR, int skipFirst_n)
{
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);

    //Market Cap
    Matrix MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Vector Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);

    cout << "Files was Loaded.\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Period");
    folderName = "Data/Output/Period/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());

    string dirName = folderName + "/Period";
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor"};
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
}

void Era_Period_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)   //TODO: Make this
{
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //Market Cap
    Matrix MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Vector Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);

    cout << "Era_PrePost_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Era_Period_PrePost");
    folderName = "Data/Output/Era_Period_PrePost/" + folderName;

    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor"};

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
            twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr - 1], iPeriods[periodNr], {Inflation_Factor[periodNr-1], Inflation_Factor[periodNr]});

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
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //Market Cap
    Matrix MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Vector Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);

    cout << "Era_PrePost_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Era_PrePost_Period");
    folderName = "Data/Output/Era_PrePost_Period/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/MarketCap.txt", "/infl_factor"};

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
            twoPeriod_Data = TwoPeriod_Calc(DR, MC, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr], {Inflation_Factor[periodNr-1], Inflation_Factor[periodNr]});

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