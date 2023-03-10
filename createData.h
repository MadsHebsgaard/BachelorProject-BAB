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


#pragma once


void DuoPeriod_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
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

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_Eras, true);

    cout << "DuoPeriod_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Double");
    folderName = "Data/Output/Double/" + folderName;
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;
    vector<vector<Matrix>> tP_DataSet(Era_List.size(),vector<Matrix>(2,Matrix(0,Vector(0))));

    for (int Era = 0; Era < Era_List.size(); ++Era)
    {
        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 0; Period < Era_List[Era].size()-1; Period++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, Era_List[Era][Period], Era_List[Era][Period + 1]);
            push_back(tP_DataSet[Era], twoPeriod_Data);
        }
        //Create directory and path for era
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());
        string preDirName = dirName + "/Pre_Period";
        string periodDirName = dirName + "/Period";

        //Create files in directory for the era
        vector<string> paths {preDirName, periodDirName};
        vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
        for (int i = 0; i < paths.size(); ++i) {
            mkdir(paths[i].c_str());
            for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(paths[i] + fileNames[file],tP_DataSet[Era][i][file]);
        }
    }
}
void SingleEra_Calculations(string folderName, double max_ratio, int minTradingDays, int n_Eras, Matrix DR)
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
    mkdir("Data/Output/Single");
    folderName = "Data/Output/Single/" + folderName;
    string dirName = folderName + "/Era";
    mkdir(folderName.c_str());
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

    for (int i = 0; i < Era_List.size(); ++i)
    {
        //Vector beta(0), alpha(0), PERMNO(0), akk_return(0), akk_sp500(0), akk_riskFree(0);    //TODO: make it vectors

        for (auto & iPeriod : Era_List[i])
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, 130);
            push_back(DataSet[i], beta_alpha_return);
        }
        // Creating a directory for the Era and each file in Era
        string periodDirName = dirName + "_" + to_string(i + 1);
        mkdir(periodDirName.c_str());
        for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(periodDirName + fileNames[file], DataSet[i][file]);
    }
}
void SinglePeriod_Calculations(string folderName, double max_ratio, int minTradingDays, Matrix DR, int skipFirst_n)
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