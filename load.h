#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "calculations.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

void HowToGetStarted()
{
    mkdir("Data");
    mkdir("Data/Input");
    string exoDirName = "Data/Input/Exo_Files";
    mkdir(exoDirName.c_str());
    cout << "\nYou need to put the following .txt files in the folder \"" << exoDirName << "\":" << endl;
    cout << "   - " << "DR\n";
    cout << "   - " << "DailyYearlyRiskFreeReturn\n";
    cout << "   - " << "sp500\n";
    cout << "   - " << "DateList\n";
    cout << "   - " << "Mth_PrevCap\n";

    cout << "When all the files are present, Run \"Process_Files()\"." << endl << endl;
}

Vector Load_Vector(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Vector(0);   }
    Vector v;
    double e;
    while(!fil.eof())
    {
        fil >> e;
        v.push_back(e);
    }
    return v;
}
Matrix Load_Matrix(const string& fn)
{
    //Forbind til fil
    ifstream fil(fn);
    if(!fil)    { cout << "Load_Matrix: Åbning af filen " << fn << " mislykkedes.\n" ; return Matrix(0);   }

    //Lav matrix struktur
    size_t m,n;
    fil >> m;
    fil >> n;
    Matrix A(m);
    for (auto& r: A)    r.resize(n);

    //Indhent matrix elementer
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)  fil >> A[i][j];
    }
    return A;
}
Intor Load_Intor(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_Intor: Could not read the file " << fn << ".\n";  return Intor(0);   }
    Intor v;
    int e;
    while(!fil.eof())
    {
        fil >> e;
        v.push_back(e);
    }
    //cout << "\nLoad_Intor: Sucssesfully loaded " << v.size() << " elements from " << fn << " to Intor.\n";
    return v;
}
Intrix Load_Dates(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }

    Intrix A;
    string junk_string;
    int IDnew, junk_number, i=0;
    vector<int> Stock(3);
    fil >> junk_string;    fil >> junk_string;    fil >> junk_string;

    while(fil >> IDnew)
    {
        Stock[0] = IDnew;
        fil >> Stock[1];
        fil >> Stock[2];

        while(IDnew == Stock[0] && !fil.eof())
        {
            fil >> IDnew;
            fil >> junk_number;
            fil >> junk_number;
        }
        //cout << "ID = " << Stock[0] << ", i = " << i << endl;
        //if(i % 500 == 0) cout << "Stock[0] = " << Stock[0] << ", i = " << i << endl;
        A.push_back(Stock);
        i++;
    }
    //cout << "All " << i << " dates was loaded sucsessfully.";
    return A;
}
Matrix Load_DR(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Matrix(0);   }
    Matrix DR;
    Vector Stock;
    int IDnew, date, date_old, i=0;
    double r, r_;
    string junk;
    int Blanc_r = -2;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;

    cout << "Loaded 0/~36150 = 0%.\n";

    while(!fil.eof())
    {
        Stock = Vector(0);
        Stock.push_back(IDnew);

        {
            if(fil.eof()) break;
            fil >> date;
            fil >> r;
            if(r < 10000)   fil >> IDnew;
            else
            {
                IDnew = r;
                r = Blanc_r;     //When blanc/no return
            }
            Stock.push_back(r);
        }
        while(true)
        {
            if(IDnew != Stock[0]) break;
            else if(fil.eof()) break;
            date_old = date;
            fil >> date;

            while(date_old == date)
            {
                fil >> r;
                if(r < 10000)  fil >> IDnew;
                else    {IDnew = r; r = Blanc_r;}
                if(IDnew != Stock[0])    break;
                date_old = date;
                fil >> date;
                if(fil.eof()) break;
            }
            if(IDnew != Stock[0] or fil.eof())    break;
            fil >> r_;
            if(r_ < 10000)
            {
                r = r_;
                fil >> IDnew;
            }
            else
            {
                r = Blanc_r;
                IDnew = r_;
            }
            Stock.push_back(r);
        }
        DR.push_back(Stock);
        i++;
        if(i % 500 == 0) cout << "Loaded " << i << "/~36150 = " << i/361.5 << "%.\n";
        if(i==max)   break;
    }
    cout << "Loaded " << i << "/~36150 = " << i / 361.5 << "%.\n";
    return DR;
}
Intrix Load_Dates_from_DR(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_Dates_from_DR: Could not read the file " << fn << ".";  return Intrix(0);   }
    Intrix Dates;
    Intor Stock;
    int IDnew, date, i=0;
    double r;
    string junk;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;
    cout << "Load_Dates_from_DR: Loaded 0/~37000 = 0%.\n";

    while(!fil.eof())
    {
        Stock = Intor(0);
        fil >> date;
        Stock.push_back(IDnew);
        Stock.push_back(date);
        while(!fil.eof())
        {
            fil >> r;
            if(r < 10000)   fil >> IDnew;
            else            IDnew = r;
            if(IDnew != Stock[0])   break;
            fil >> date;
        }
        Stock.push_back(date);
        Dates.push_back(Stock);
        i++;
        if(i % 500 == 0) cout << "Load_Dates_from_DR: Loaded " << i << "/~36150 = " << i/370 << "%.\n";
    }
    cout << "Load_Dates_from_DR: Loaded " << i << "/~37000 = " << i / 370 << "%.\n";
    return Dates;
}
Matrix Load_DR_Compressed(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_DR_Compressed: Could not read the file " << fn << ".";  return Matrix(0);   }

    Vector Stock;
    Matrix DR;

    DR.reserve(40000);
    int n, i=0;

    while(fil >> n)
    {
        Stock = Vector(n);
        for (int j = 0; j < n; ++j)
        {
            fil >> Stock[j];
        }
        DR.push_back(Stock);
        if(i % 500 == 0)    cout << "ID = " << Stock[0] << " i = " << i << endl;
        i++;
        if (i == max) break;
    }
    cout << "Load_DR_Compressed: Sucssesfully loaded " << i << " stocks to DR from " << fn << ".\n";
    return DR;
}
Matrix Load_MC_Compressed(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_DR_Compressed: Could not read the file " << fn << ".";  return Matrix(0);   }
    Matrix MC;

    int n;
    while(fil >> n)
    {
        Vector stock(n);
        for (int j = 0; j < n; ++j)
        {
            fil >> stock[j];
        }
        MC.push_back(stock);
    }
    return MC;
}
Intrix Load_Intrix(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_Intrix: Could not read the file " << fn << ".\n";  return Intrix(0);   }

    if(max < 1) max = 999999;

    int m, n;
    fil >> m;
    fil >> n;
    Intrix A(min(m, max), Intor(n));

    for (int i = 0; i < min(m, max); ++i)
        for (int j = 0; j < n; ++j)
            fil >> A[i][j];

    /*if ( A[0].size() == A[1].size() || A[0].size() == A[A.size()-1].size() )
        cout << "Load_Intrix: Sucessfully loaded (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
    else    cout << "Load_Intrix: Sucssesfully loaded (" << A.size() << " x X) to Intrix.\n";
*/
    return A;
}
Intrix Load_StockDays_from_DR(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }
    Intrix Stockdays;
    Intor days;
    int IDnew, date, date_old, i=0;
    double r_;
    string junk;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;
    cout << "Loaded 0/36148 = 0%.\n";

    while(!fil.eof())
    {
        days = Intor(0);
        days.push_back(IDnew);

        fil >> date;
        days.push_back(date);
        fil >> r_;
        if(r_<10000)  fil >> IDnew;
        else          IDnew = r_;

        while(true)
        {
            if(IDnew != days[0]) break;
            else if(fil.eof()) break;
            date_old = date;
            fil >> date;

            while(date_old == date)
            {
                fil >> r_;
                if(r_ < 10000)  fil >> IDnew;
                else    IDnew = r_;
                if(IDnew != days[0])    break;
                date_old = date;
                fil >> date;
                if(fil.eof()) return Stockdays;
            }
            if(IDnew != days[0])    break;
            days.push_back(date);
            fil >> r_;
            if(r_ < 10000)
            {
                fil >> IDnew;
            }
            else
            {
                IDnew = r_;
            }
        }
        Stockdays.push_back(days);
        i++;
        if(i % 500 == 0) cout << "Loaded " << i << "/36148 = " << i/361.48 << "%.\n";
        if(i == max)   break;
    }
    cout << "Loaded " << i << "/36148 = " << i / 361.48 << "%.\n";
    return Stockdays;
}
Intrix Load_StockDays_Compressed(const string& fn, bool With_Id, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }
    cout << "Loading " << max << " StockDays from compressed file\n";
    Intor days;
    Intrix StockDays;
    int n, i=0;

    while(fil >> n)
    {
        if(!With_Id)     n--;
        days = Intor(n);
        int j = 0;
        if(With_Id) {fil >> days[0];    j++;}
        while (j < n)
        {
            fil >> days[j];
            j++;
        }
        StockDays.push_back(days);
        if(i % 2000 == 0)    cout << "Loaded up to ID: " << days[0] << ", i = " << i << endl;
        i++;
        if (i == max) break;
    }
    cout << "Sucssesfully loaded " << i << " stocks to StockDays.\n";
    return StockDays;
}
Matrix Load_Mth_MarketCap(string fn, int Factor)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".\n";  return Matrix(0);   }
    //fil >> std::setprecision(2);

    Matrix MC(0);
    string junk;
    int IDnew;
    int IDold;
    int date;
    double PrevCap;
    string string_number;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;
    fil >> date;

    IDold = IDnew;
    Vector Stock(2);
    Stock[0] = IDnew;

    bool ZeroMonth = true;
    bool LaterMonth = false;

    while (true)
    {
        if (IDnew != IDold) {
            MC.push_back(Stock);
            Stock = Vector(2);
            Stock[0] = IDnew;
            Stock[1] = -2;   //nan("") is maybe better
            IDold = IDnew;
            //PrevCap = -0.002;
            ZeroMonth = true;
            LaterMonth = false;
        }
        fil >> string_number;

        if (string_number.contains('.')) {
            ZeroMonth = false;
            PrevCap = stod(string_number); //string_number = MthPrevCap
            if(fil.eof()){  Stock.push_back(PrevCap); MC.push_back(Stock); }    //Rare case if eof
            fil >> IDnew;
        }
        else    IDnew = stoi(string_number);    //string_number = IDnew;
        if(!ZeroMonth)  {
            if(!LaterMonth)
                Stock[1] = date/100;
            Stock.push_back(PrevCap * 1000/Factor);
            LaterMonth = true;
        }
        fil >> date;
        if(fil.eof())   {MC.push_back(Stock); break;}
    }
    cout << "Load_Monthly_MarketCap is done...\n";
    return MC;
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
    logMessage.push_back("max_ratio = "+to_string(max_ratio));
    logMessage.push_back("minTradingDays = "+to_string(minTradingDays));
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //DR with condition for inclusion
    DR = Edit_DR(DR, max_ratio, minTradingDays);
    logMessage.push_back("DR.size() = "+to_string(DR.size()));

    //Market Cap
    MC = Load_MC_Compressed("Data/Input/Processed_Files/MarketCap_yr.txt");
    MC = Remove_Missing_ID(MC, Matrix_Column(DR, 0));

    //Inflation factor
    Inflation_Factor = Load_Vector("Data/Input/Processed_Files/Inflation_Factor.txt");

    //iDates with same stocks as DR
    iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID(iDates, Matrix_Column(DR, 0));

    //Load sp500 and riskFree returns & SP500 dates
    sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);

    //Create files and folderName = folderPath
    mkdir("Data/Output");
    cout << methodName << "_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir(("Data/Output/" + methodName).c_str());
    folderName = "Data/Output/" + methodName + "/" + folderName;
    std::filesystem::remove_all(folderName.c_str());
    mkdir(folderName.c_str());
}