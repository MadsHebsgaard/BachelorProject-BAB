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

void Save(const string& fn, const Vector& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Save: Åbning af filen " << fn << " mislykkedes.\n" ; return;   }
    fil << fixed << setprecision(8);
    for(double e : v)   {
        if(floor(e) == e)   fil << setprecision(0);
        else                fil << setprecision(8);
        fil << e << endl;
    }
}
void Save(const string& fn, const Intor& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(int e : v)   fil << e << endl;
}
void Save(const string& fn, Intrix A)
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
    //cout << "Save: Sucessfully saved (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
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
void Compress_MC(const string& fn, Matrix MC, int factor)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    fil << fixed << setprecision(0+floor(log10(factor)));

    int i=0;
    for (auto & stock : MC)
    {
        fil << stock.size() << " " << stock[0] << " " << stock[1];
        for (int j = 2; j < stock.size(); ++j)  {
            if((j-2) % 10 == 0) fil << endl;
            fil << stock[j] << " ";
        }
        fil << endl << endl;
    }
}
void Process_Files() {
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";
    const std::vector<std::string> filenames = {"sp500.txt", "DR.txt", "DateList.txt", "DailyYearlyRiskFreeReturn.txt, Mth_PrevCap.txt"};

    if (areFilesExistInDirectory(filenames, Exo_FilePath))
    {
        cout << "\nCreating all necessary files in " << Proccessed_FilePath << ":\n\n";
        mkdir("Data/Input/Processed_Files");
        int max = 99999999;

        //Daily Return on each stock compressed
        Matrix DR = Load_DR(Exo_FilePath + "DR.txt", max);  //TODO: uncomment
        Compress_DR(Proccessed_FilePath + "DR_Compressed.txt", DR);   //TODO: uncomment
        cout << "Created " << Proccessed_FilePath << "DR_Compressed.txt\n";

        //Each Stock's lifespan
        Intrix DR_Dates = Load_Dates_from_DR(Exo_FilePath + "DR.txt");
        Save(Proccessed_FilePath + "DR_Dates.txt", DR_Dates);
        cout << "Created " << Proccessed_FilePath << "DR_Dates.txt\n";
        //Intrix DR_Dates = Load_Intrix(Proccessed_FilePath+"DR_Dates.txt",-1); //if DR_Dates already exists

        //DR_iDates, Stock's lifespan in index values starting from (0) the first recorded data date.
        Intor DateList = Load_Intor(Exo_FilePath + "DateList.txt");
        Intrix DR_iDates = Dates_to_iDates(DR_Dates, DateList, 1);
        Save(Proccessed_FilePath + "DR_iDates.txt", DR_iDates);
        cout << "Created " << Proccessed_FilePath << "DR_iDates.txt\n";

        //iPeriods
        Intrix iPeriods = Yearly_iPeriods(DateList);
        Save(Proccessed_FilePath + "iPeriods.txt", iPeriods);
        cout << "Created " << Proccessed_FilePath << "iPeriods.txt\n";

        //DailyDailyRFR
        Vector DailyYearlyRFR = Load_Vector(Exo_FilePath + "DailyYearlyRiskFreeReturn.txt");
        Vector DailyDailyRFR = DailyYearly_to_DailyDaily_Return(DailyYearlyRFR, iPeriods);
        Save(Proccessed_FilePath + "riskFreeReturn.txt", DailyDailyRFR);
        cout << "Created " << Proccessed_FilePath << "riskFreeReturn.txt\n";

        //Market Cap
        int factor = 1; //todo: needs to be 1 atm (Compress_MC)
        Matrix MarketCap_Mth = Load_Mth_MarketCap("Data/Input/Exo_files/Mth_PrevCap.txt", factor);
        Matrix MarketCap_yearly = MarketCap_Monthly_to_Yearly(MarketCap_Mth);
        Compress_MC("Data/Input/Processed_Files/MarketCap_yr.txt", MarketCap_yearly, factor);

        //Other files, usefull for testing
        //Intrix StockDays = Load_StockDays_from_DR("DR.txt", max);
        //Compress_DR_StockDays("DR_StockDays.txt", StockDays);    } else {
    }
    else
        HowToGetStarted();
}
void LogFile(string filePath, vector<string> Message)
{
    ofstream fil(filePath + "/LogFile.txt");
    for(auto& line:Message)
        fil << line << endl;
}