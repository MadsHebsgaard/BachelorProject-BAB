#include <vector>
#include <iostream>
#include <fstream>
#include "calculations.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once

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
    fil.close();
}
void Save(const string& fn, const Intor& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(int e : v)   fil << e << endl;
    fil.close();
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
    fil.close();
}
void Save(const string& fn, Matrix A)
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
    fil.close();
}
void Compress_DR(const string& fn, Matrix DR)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    int i=0;
    for (auto & stock : DR)
    {
        int stock_size = lastNumberOverK(stock, -1.5);
        if (stock_size < 2) continue;   //Skip stock if all returns are missing '-2'
        fil << endl << endl << stock_size << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock_size; ++j)  fil << " " << stock[j];

        if(i % 500 == 0) cout << "Compressed " << i << "/" << DR.size() << " = " << (100*i)/DR.size() << "%.\n";    i++;
    }
    cout << "Compress_DR: Sucssesfully Compressed/Saved dily 'returns' of " << DR.size() << "  stocks to " << fn << ".\n";
    fil.close();
}
void Save_RS_Dates(const string& fn, Matrix DR)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    int i=0;
    for (auto & stock : DR)
    {
        int stock_size = lastNumberOverK(stock, -1.5);
        if (stock_size < 2) continue;   //Skip stock if all returns are missing '-2'
        fil << endl << endl << stock_size << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock_size; ++j)  fil << " " << stock[j];

        if(i % 500 == 0) cout << "Compressed " << i << "/" << DR.size() << " = " << (100*i)/DR.size() << "%.\n";    i++;
    }
    cout << "Compress_DR: Sucssesfully Compressed/Saved dily 'returns' of " << DR.size() << "  stocks to " << fn << ".\n";
    fil.close();
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
    fil.close();
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
    fil.close();
}

void saveLogFile(string filePath, vector<string> Message)
{
    ofstream fil(filePath + "/saveLogFile.txt");
    for(auto& line:Message)
        fil << line << endl;
    fil.close();
}
void save_logfile(auto start, vector<string> logMessage, string folderName)
{
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    logMessage.push_back("finished computation at " + string(ctime(&end_time)));
    logMessage.push_back("elapsed time: " + to_string(elapsed_seconds.count()) + "s");
    saveLogFile(folderName, logMessage);    //Add more information to logMessage
    cout << "elapsed time: " + to_string(elapsed_seconds.count()) + "s\n";
}
void Save_CSV(string fn, Matrix A, vector<string> headers)
{
    //todo: Make it so that only if '/' and/or ".txt" is present, they are removed
    remove_char(headers, '/');
    remove_substring(headers, ".txt");

    char delim = ';';
    ofstream fil(fn);
    int n = headers.size();
    for (int i = 0; i<n-1; ++i) {
        fil << headers[i] << delim;
    }
    fil << headers[n-1] << endl << fixed;

    for (int j = 0; j<A[0].size(); ++j){
        for (int i = 0; i<A.size()-1; ++i)
        {
            A[i][j] == round(A[i][j]) ? fil << setprecision(0) : fil << setprecision(8);
            fil << A[i][j] << delim;
        }
        A[A.size()-1][j] == round(A[A.size()-1][j]) ? fil << setprecision(0) : fil << setprecision(8);
        fil << A[A.size()-1][j] << endl;
    }
    fil.close();
}
void Save_TwoDataSet_CSV(string fn, vector<Matrix> T3, vector<string> headers, vector<string> header_nr)
{
    //todo: Make it so that only if '/' and/or ".txt" is present, they are removed
    remove_char(headers, '/');
    remove_substring(headers, ".txt");

    char delim = ';';
    ofstream fil(fn);
    int max_variable = T3[0].size();
    int max_observation = T3[0][0].size();
    int s_max = T3.size();

    for(int s=0; s<s_max; s++)
        for (int i = 0; i<max_variable; ++i)
            if(s+1<s_max || i<max_variable-1)
                fil << headers[i] << header_nr[s] << delim;
    fil << headers[max_variable-1] << header_nr[s_max-1] << endl << fixed;

    for (int i = 0; i<max_observation; ++i)
    {
        for(int s = 0; s<s_max; s++)
            for (int j = 0; j<max_variable; ++j)
            {
                T3[s][j][i] == round(T3[s][j][i]) ? fil << setprecision(0) : fil << setprecision(8);
                if(s<s_max-1 || j<max_variable-1)
                        fil << T3[s][j][i] << delim;
                else    fil << T3[s_max-1][j][i] << endl;
            }
    }
    fil.close();
}
