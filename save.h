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

#include "calculations.h"

void Save_Vector(const string& fn, const Vector& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Save_Vector: Åbning af filen " << fn << " mislykkedes.\n" ; return;   }
    cout << fixed << setprecision(15);
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
    //cout << "Save_Intrix: Sucessfully saved (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
}
void Compress_DR(const string& fn, Matrix DR)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    int i=0;
    for (auto & stock : DR)
    {
        fil << endl << endl << stock.size() << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock.size(); ++j)  fil << " " << stock[j];

        if(i % 500 == 0) cout << "Compressed " << i << "/" << DR.size() << " = " << (100*i)/DR.size() << "%.\n";    i++;
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

void Process_Files() {
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";
    const std::vector<std::string> filenames = {"sp500.txt", "DR_No_Ticker.txt", "DateList.txt", "DailyYearlyRiskFreeReturn.txt"};

    if (areFilesExistInDirectory(filenames, Exo_FilePath))
    {
        cout << "\nCreating all necessary files in " << Proccessed_FilePath << ":\n\n";
        mkdir("Data/Input/Processed_Files");
        int max = 99999999;

        //Daily Return on each stock compressed
        Matrix DR = Load_DR(Exo_FilePath + "DR_No_Ticker.txt", max);  //TODO: uncomment
        Compress_DR(Proccessed_FilePath + "DR_Compressed.txt", DR);   //TODO: uncomment
        cout << "Created " << Proccessed_FilePath << "DR_Compressed.txt\n";

        //Each Stock's lifespan
        Intrix DR_Dates = Load_Dates_from_DR(Exo_FilePath + "DR_No_Ticker.txt");
        Save_Intrix(Proccessed_FilePath + "DR_Dates.txt", DR_Dates);
        cout << "Created " << Proccessed_FilePath << "DR_Dates.txt\n";
        //Intrix DR_Dates = Load_Intrix(Proccessed_FilePath+"DR_Dates.txt",-1); //if DR_Dates already exists

        //DR_iDates, Stock's lifespan in index values starting from (0) the first recorded data date.
        Intor DateList = Load_Intor(Exo_FilePath + "DateList.txt");
        Intrix DR_iDates = Dates_to_iDates(DR_Dates, DateList, 1);
        Save_Intrix(Proccessed_FilePath + "DR_iDates.txt", DR_iDates);
        cout << "Created " << Proccessed_FilePath << "DR_iDates.txt\n";

        //iPeriods
        Intrix iPeriods = Yearly_iPeriods(DateList);
        Save_Intrix(Proccessed_FilePath + "iPeriods.txt", iPeriods);
        cout << "Created " << Proccessed_FilePath << "iPeriods.txt\n";

        //DailyDailyRFR
        Vector DailyYearlyRFR = Load_Vector(Exo_FilePath + "DailyYearlyRiskFreeReturn.txt");
        Vector DailyDailyRFR = DailyYearly_to_DailyDaily_Return(DailyYearlyRFR, iPeriods);
        Save_Vector(Proccessed_FilePath + "riskFreeReturn.txt", DailyDailyRFR);
        cout << "Created " << Proccessed_FilePath << "riskFreeReturn.txt\n";

        //Other files, usefull for testing
        //Intrix StockDays = Load_StockDays_from_DR("DR_No_Ticker.txt", max);
        //Compress_DR_StockDays("DR_StockDays.txt", StockDays);    } else {
    }
    else
        HowToGetStarted();
}