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
//DR_No_Ticker

void Era_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_Eras, true);

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

    for (int Era = 0; Era < Era_List.size(); ++Era)
    {
        vector<Matrix> tP_DataSet(vector<Matrix>(2,Matrix(0,Vector(0))));

        if(Era < rest)  eraLength = minEraLength+1;
        else            eraLength = minEraLength;

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 1; Period <= eraLength; Period++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr]);
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
        vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
        for (int prePost = 0; prePost < paths.size(); ++prePost) {
            mkdir(paths[prePost].c_str());
            for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(paths[prePost] + fileNames[file], tP_DataSet[prePost][file]);
        }
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
void Era_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
{
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    cout << "Files was Loaded.\n\n";

    Matrix beta(n_Eras, Vector(0)), alpha(n_Eras, Vector(0)), PERMNO(n_Eras, Vector(0));
    Matrix akk_return(n_Eras, Vector(0)), akk_sp500(n_Eras, Vector(0)), akk_riskFree(n_Eras, Vector(0));
    Matrix beta_alpha_return;

    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_Eras, true);

    vector<Matrix> DataSet(Era_List.size());

    mkdir("Data/Output");
    mkdir("Data/Output/Era");
    folderName = "Data/Output/Era/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());

    string dirName = folderName + "/Era";
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

    for (int i = 0; i < Era_List.size(); ++i)
    {
        //Vector beta(0), alpha(0), PERMNO(0), akk_return(0), akk_sp500(0), akk_riskFree(0);    //TODO: make it vectors

        for (auto & iPeriod : Era_List[i])
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, minTradingDays);
            push_back(DataSet[i], beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(i + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(periodDirName + fileNames[file], DataSet[i][file]);
    }
}
void Period_Calculations(string folderName, double max_ratio, int minTradingDays, Matrix DR, int skipFirst_n)
{
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR, max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

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
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
    Matrix Data;

    for (int period = skipFirst_n; period < iPeriods.size(); ++period)
    {
        //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
        Data = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriods[period], 130);

        // Creating a directory and files for the period
        string periodDirName = dirName + "_" + to_string(period + 1 - skipFirst_n);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(periodDirName + fileNames[file], Data[file]);
    }
}

void Era_Period_PrePost_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)   //TODO: Make this
{
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";


    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);
    logMessage.push_back("DR.size() = "+to_string(DR.size()));
    logMessage.push_back("DR_ny.size() = "+to_string(DR_ny.size()));

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

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
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

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
            if(eraPeriod==61)
            {
                cout << "her";
            }
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, iPeriods[periodNr -1], iPeriods[periodNr]);

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
                for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(prePost[i] + fileNames[file], twoPeriod_Data[i][file]);
            }
            periodNr++;
        }
        cout << "Era " << Era << " is done.\n";
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
void Era_PrePost_Period_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)   //TODO: Make this
{
    vector<string> logMessage = {"max_ratio = "+to_string(max_ratio),"minTradingDays = "+to_string(minTradingDays)};
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

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
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

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
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, iPeriods[periodNr-1], iPeriods[periodNr]);

            //Create dirs and then files
            string periodPre_DirName = preDirName + "/Period_" + to_string(eraPeriod);
            string periodPost_DirName = postDirName + "/Period_" + to_string(eraPeriod);

            mkdir(periodPre_DirName.c_str());
            mkdir(periodPost_DirName.c_str());

            vector<string> prePost {periodPre_DirName, periodPost_DirName};
            for (int i = 0; i < prePost.size(); ++i) {
                mkdir(prePost[i].c_str());
                for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(prePost[i] + fileNames[file], twoPeriod_Data[i][file]);
            }
            periodNr++;
        }
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}

void BackTesting(Tensor5 Data, int eraNr)
{


    //vector<Matrix> Low_Portfolio(Data[eraNr][0].size(), Matrix(stockInf, Vector(0)));
    //vector<Matrix> High_Portfolio(Data[eraNr][0].size(), Matrix(stockInf, Vector(0)));

    //Creating the portfolio with all the stocks
    //for (int eraNr = 0; eraNr < Data.size(); ++eraNr) {
    /*for (int periodNr = 0; periodNr < Data[eraNr][0].size(); ++periodNr) {
        for (int stockNr = 0; stockNr < Data[eraNr][0][periodNr][0].size(); ++stockNr) {

            //Low beta portefolio
            if(Data[eraNr][0][periodNr][0][stockNr] < betaLow)
            {
                for (int inform = 0; inform < stockInf; ++inform)
                    Low_Portfolio[periodNr][inform].push_back(Data[eraNr][1][periodNr][inform][stockNr]);
            }
            //High beta portefolio
            else if(Data[eraNr][0][periodNr][0][stockNr] > betaLow)
            {
                for (int inform = 0; inform < stockInf; ++inform)
                    High_Portfolio[periodNr][inform].push_back(Data[eraNr][1][periodNr][inform][stockNr]);
            }
        }
    }*/
    //}

    double betaLow = 0.5;
    double betaHigh = 1.25;
    int stockInf = Data[0][0][0].size(); //Permno is excluded
    eraNr--;

    int iBeta = 0;
    int iReturn = 1;
    int iAlpha = 2;
    int iSP500 = 3;
    int iRiskFree = 4;
    int iCount = 5;
    stockInf = Data[0][0][0].size()-1; //Permno is excluded

    int periodCount = Data[eraNr].size();

    Matrix lowPortfolio(periodCount, Vector(stockInf, 0));
    Matrix highPortfolio(periodCount, Vector(stockInf, 0));
    Matrix BAB(periodCount, Vector(stockInf, 0));

    Matrix preBeta(periodCount, Vector(2,0));
    Matrix Count(periodCount, Vector(2,0));

    Vector BAB_Average(stockInf+1, 0);

    //Creating the portfolio as if it was a stock (ETF)
    for (int periodNr = 0; periodNr < periodCount; ++periodNr) {
        for (int stockNr = 0; stockNr < Data[eraNr][periodNr][0][0].size(); ++stockNr) {

            //Low beta portefolio
            if(Data[eraNr][periodNr][0][0][stockNr] < betaLow)
            {
                preBeta[periodNr][0] += Data[eraNr][periodNr][0][0][stockNr];
                Count[periodNr][0]++;
                for (int inform = 0; inform < stockInf; ++inform) {
                    lowPortfolio[periodNr][inform] += Data[eraNr][periodNr][1][inform][stockNr];
                }
            }
                //High beta portefolio
            else if(Data[eraNr][periodNr][0][0][stockNr] > betaLow)
            {
                preBeta[periodNr][1] += Data[eraNr][periodNr][0][0][stockNr];
                Count[periodNr][1]++;
                for (int inform = 0; inform < stockInf; ++inform) {
                    highPortfolio[periodNr][inform] += Data[eraNr][periodNr][1][inform][stockNr];
                }
            }
        }
        preBeta[periodNr][0] /= Count[periodNr][0];
        preBeta[periodNr][1] /= Count[periodNr][1];

        for (int info = 0; info < stockInf; ++info) {
            lowPortfolio[periodNr][info] /= Count[periodNr][0];     //Average
            highPortfolio[periodNr][info] /= Count[periodNr][1];    //Average
        }
        double LowAmmount = preBeta[periodNr][1] / (preBeta[periodNr][0]+preBeta[periodNr][1]);
        for (int info = 0; info < stockInf; ++info) {
            BAB[periodNr][info] = LowAmmount * lowPortfolio[periodNr][info] - (1-LowAmmount) * highPortfolio[periodNr][info];
            BAB_Average[info] += BAB[periodNr][info];
        }
        BAB_Average[iCount] += Count[periodNr][0] + Count[periodNr][1];
    }
    for (int info = 0; info < stockInf+1; ++info)   BAB_Average[info] /= periodCount;
    return;
}