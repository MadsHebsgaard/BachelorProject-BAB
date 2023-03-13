#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once
#include "load.h"

int Find_date_integer(int date, Intor& Date_list);


void day_est_checker(int max)
{
    Intrix Dates = Load_Intrix("Dates_DR_Compressed.txt", max);
    Intor SP500_Dates = Load_Intor("List_of_dates.txt");

    if(max < 0 or max > Dates.size()) max = Dates.size();

    int True_date, Est_date;
    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][1];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][1], SP500_Dates)];
        if(True_date < Est_date)    cout << True_date << " <1 " << Est_date << ",  i = " << i << endl;
    }

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][2];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][2], SP500_Dates)];
        if(True_date < Est_date)    cout << True_date << " <2 " << Est_date << ",  i = " << i << endl;
    }

    cout << " --------------------------------------------- " << endl << endl;

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][1];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][1], SP500_Dates)];
        if(True_date > Est_date)    cout << True_date << " >1 " << Est_date << ",  i = " << i << endl;
    }

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][2];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][2], SP500_Dates)];
        if(True_date > Est_date)    cout << True_date << " >2 " << Est_date << ",  i = " << i << endl;
    }
}
void Date_Checker(int max)
{
    if (max < 0)    max = 36146;
    Intrix StockDays = Load_StockDays_Compressed("DR_StockDays.txt",true, max);
    Intor Dates = Load_Intor("SP_Dates.txt");
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);

    //78094 & 15451

    cout << endl;

    for (int i = 0; i < max; ++i)
    {
        //int x = 300;   //1000  Last is unequal 19961105 19961106 78094   20070625
        for(int x = 0; x < iDates[i][2]-iDates[i][1]; ++x)
        {
            if(iDates[i][2] - iDates[i][1] >= x)
            {
                if(Dates[iDates[i][1] + x] != StockDays[i][1 + x])
                {
                    cout << "Last is unequal " << Dates[iDates[i][1] + x];
                    cout << " " << StockDays[i][1 + x] << " " << StockDays[i][0] << "   ";
                    cout << Dates[iDates[i][2]] << endl;
                }
            }
            else
            {
                //cout << setw(5) << iDates[i][2] - iDates[i][1] + 1 << " < " << x+1 << ". ID = " << iDates[i][0] << endl;
            }
        }
    }
}
void OldTestCalculations(int max, double max_ratio, int minTradingDays)
{
    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,max_ratio, minTradingDays);

    Intrix iPeriods = Load_Intrix("Yearly_iPeriods_1926_1949.txt", -1);
    //Intrix iPeriods = Load_Intrix("Yearly_iPeriods_1950_1973.txt", -1);
    //Intrix iPeriods = Load_Intrix("Yearly_iPeriods_1974_1997.txt", -1);
    //Intrix iPeriods = Load_Intrix("Yearly_iPeriods_1998_2021.txt", -1);
    //Intrix iPeriods = Load_Intrix("Yearly_iPeriods.txt", -1);

    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    cout << "Files was Loaded.\n\n";

    Vector beta(0), alpha(0), akk_return(0), PERMNO(0), akk_sp500(0), akk_riskFree(0);
    Matrix beta_alpha_return;
    for (auto & iPeriod : iPeriods)
    {
        //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
        beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, 130);

        //Store Data
        beta.insert(beta.end(), beta_alpha_return[0].begin(), beta_alpha_return[0].end());
        alpha.insert(alpha.end(), beta_alpha_return[1].begin(), beta_alpha_return[1].end());
        akk_return.insert(akk_return.end(), beta_alpha_return[2].begin(), beta_alpha_return[2].end());
        PERMNO.insert(PERMNO.end(), beta_alpha_return[3].begin(), beta_alpha_return[3].end());
        akk_sp500.insert(akk_sp500.end(), beta_alpha_return[4].begin(), beta_alpha_return[4].end());
        akk_riskFree.insert(akk_riskFree.end(), beta_alpha_return[5].begin(), beta_alpha_return[5].end());
    }
    //Save Data
    Save_Vector("Data_Test/beta.txt",beta);
    Save_Vector("Data_Test/alpha.txt",alpha);
    Save_Vector("Data_Test/akk_return.txt",akk_return);
    Save_Vector("Data_Test/PERMNO.txt",PERMNO);
    Save_Vector("Data_Test/akk_sp500.txt", akk_sp500);
    Save_Vector("Data_Test/akk_riskFree.txt", akk_riskFree);
}

void TestCalculations(string folderName, int max, double max_ratio, int minTradingDays, int n_periods)
{
    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");

    cout << "Files was Loaded.\n\n";

    Matrix beta(n_periods,Vector(0)), alpha(n_periods,Vector(0)), PERMNO(n_periods,Vector(0));
    Matrix akk_return(n_periods,Vector(0)), akk_sp500(n_periods,Vector(0)), akk_riskFree(n_periods,Vector(0));
    Matrix beta_alpha_return;

    Intrix iPeriods = Load_Intrix("Yearly_iPeriods.txt", -1);
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_periods, true);

    folderName = "Output_Data/" + folderName;
    mkdir(folderName.c_str());
    for (int i = 0; i < Era_List.size(); ++i)
    {
        for (auto & iPeriod : Era_List[i])
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, 130);

            //Store Data
            beta[i].insert(beta[i].end(), beta_alpha_return[0].begin(), beta_alpha_return[0].end());
            alpha[i].insert(alpha[i].end(), beta_alpha_return[1].begin(), beta_alpha_return[1].end());
            akk_return[i].insert(akk_return[i].end(), beta_alpha_return[2].begin(), beta_alpha_return[2].end());
            PERMNO[i].insert(PERMNO[i].end(), beta_alpha_return[3].begin(), beta_alpha_return[3].end());
            akk_sp500[i].insert(akk_sp500[i].end(), beta_alpha_return[4].begin(), beta_alpha_return[4].end());
            akk_riskFree[i].insert(akk_riskFree[i].end(), beta_alpha_return[5].begin(), beta_alpha_return[5].end());
        }
        // Creating a directory for the Era_List
        string dirName = folderName + "/Era";
        string periodDirName = dirName + "_" + to_string(i + 1);
        mkdir(periodDirName.c_str());

        //Save Data
        Save_Vector(periodDirName + "/beta.txt",beta[i]);
        Save_Vector(periodDirName + "/alpha.txt",alpha[i]);
        Save_Vector(periodDirName + "/akk_return.txt",akk_return[i]);
        Save_Vector(periodDirName + "/PERMNO.txt",PERMNO[i]);
        Save_Vector(periodDirName + "/akk_sp500.txt", akk_sp500[i]);
        Save_Vector(periodDirName + "/akk_riskFree.txt", akk_riskFree[i]);
    }
}



/*
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
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_Eras, true);

    cout << "Era_PrePost_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Era_PrePost_Period");
    folderName = "Data/Output/Era_PrePost_Period/" + folderName;
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;
    vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

    for (int Era = 0; Era < n_Eras; ++Era)
    {
        //Create directory and path for era + pre(period)
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());
        string preDirName = dirName + "/Pre_Period";
        string postDirName = dirName + "/Period";
        mkdir(preDirName.c_str());
        mkdir(postDirName.c_str());

        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 1; Period < Era_List[Era].size(); Period++)
        {
            //if(Era == 0 && Period == 0) continue;
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, Era_List[Era][Period-1], Era_List[Era][Period]);

            //Create dirs and then files
            string periodPre_DirName = preDirName + "/Period_" + to_string(Period);
            string periodPost_DirName = postDirName + "/Period_" + to_string(Period);

            mkdir(periodPre_DirName.c_str());
            mkdir(periodPost_DirName.c_str());

            vector<string> prePost {periodPre_DirName, periodPost_DirName};
            for (int i = 0; i < prePost.size(); ++i) {
                mkdir(prePost[i].c_str());
                for (int file = 0; file < fileNames.size(); ++file)    Save_Vector(prePost[i] + fileNames[file], twoPeriod_Data[i][file]);
            }
        }
    }
    LogFile(folderName, logMessage);    //Add more information to logMessage
}
 */