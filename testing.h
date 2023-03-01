#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once
#define BACHELOR_TESTING_H
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
