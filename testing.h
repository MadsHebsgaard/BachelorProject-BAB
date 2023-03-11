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
