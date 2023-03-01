#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once
#define BACHELOR_SAVE_H

Intrix Dates_to_iDates(Intrix Dates, Intor True_Dates, int start_j);

void Save_Vector(const string& fn, const Vector& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(double e : v)   fil << e << endl;
}
void Save_Intor(const string& fn, const Intor& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(int e : v)   fil << e << endl;
}
void Save_Intrix(const string& fn, Intrix A)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    fil << A.size() << endl;
    fil << A[0].size() << endl << endl;

    for (int i = 0; i < A.size(); ++i)
    {
        //fil << setw(10) << A[i][0];
        for (int j = 0; j < A[0].size(); ++j)
        {
            fil << setw(10) << A[i][j];
        }
        fil << endl;
    }
    cout << "Save_Intrix: Sucessfully saved (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
}
void Compress_DR(const string& fn, Matrix DR)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    for (auto & stock : DR)
    {
        fil << endl << endl << stock.size() << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock.size(); ++j)  fil << " " << stock[j];
    }
    cout << "Compress_DR: Sucssesfully Compressed/Saved dily 'returns' of " << DR.size() << "  stocks to " << fn << ".\n";
}
void Compress_DR_StockDays(const string& fn, Intrix StockDays)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    for (auto & stock : StockDays)
    {
        fil << endl << endl << stock.size() << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock.size(); ++j)  fil << " " << stock[j];
    }
    cout << "Compress_DR_StockDays: Sucssesfully Compressed/Saved daily 'dates' of " << StockDays.size() << " stocks to " << fn << ".\n";
}
void sp500_Dates_to_monthly_dDate_periods()
{
    Intor allDates = Load_Intor("SP_Dates.txt");
    int i = 0, periodStart = allDates[0];
    Intrix MonthlyPeriods(0);

    while (allDates[i] < 99999999)
    {
        if( (allDates[i] - periodStart) > 40)
        {
            MonthlyPeriods.push_back({periodStart, allDates[i-1]});
            periodStart = allDates[i];
        }
        i++;
    }
    //Save_Intrix("Monthly_dPeriods.txt", MonthlyPeriods);
    Intrix Monthly_iPeriods = Dates_to_iDates(MonthlyPeriods,allDates,0);
    Save_Intrix("Monthly_iPeriods.txt", Monthly_iPeriods);


    i = 0;
    periodStart = allDates[0];
    Intrix YearlyPeriods(0);
    while (allDates[i] < 99999999)
    {
        if( (allDates[i] - periodStart) > 4000)
        {
            YearlyPeriods.push_back({periodStart, allDates[i-1]});
            periodStart = allDates[i];
        }
        i++;
    }
    Intrix Yearly_iPeriods = Dates_to_iDates(YearlyPeriods,allDates,0);
    Save_Intrix("Yearly_iPeriods.txt", Yearly_iPeriods);
}
void Create_Files_from_DR(int max)
{
    //Need files:
    //DR_NO_Ticker.txt
    //SP_Dates
    //Period_Dates.txt  (optional if needing periods)

    //DR
    Matrix DR = Load_DR("DR_No_Ticker.txt", max);
    Compress_DR("DR_Compressed.txt", DR);

    //StockDays - Usefull while testing
    Intrix StockDays = Load_StockDays_from_DR("DR_No_Ticker.txt", max);
    Compress_DR_StockDays("DR_StockDays.txt", StockDays);

    //DR_Dates
    Intrix DR_Dates = Load_Dates_from_DR("DR_No_Ticker.txt");
    Save_Intrix("DR_Dates.txt", DR_Dates);

    //iDates
    Intor SP_Dates = Load_Intor("SP_Dates.txt");
    Intrix iDates = Dates_to_iDates(DR_Dates, SP_Dates, 1);
    Save_Intrix("DR_iDates.txt", iDates);

    //iPeriod
    Intrix dPeriod = Load_Intrix("Period_Dates.txt", -1);
    Save_Intrix("dPeriod.txt", dPeriod);

    //Intrix dPeriod = Load_Intrix("dPeriod.txt", -1);
    Intrix iPeriod = Dates_to_iDates(dPeriod, SP_Dates, 0);
    Save_Intrix("iPeriod.txt", iPeriod);

    sp500_Dates_to_monthly_dDate_periods();
}