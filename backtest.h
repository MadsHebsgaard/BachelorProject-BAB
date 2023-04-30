#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()
#include "load_output.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;

#pragma once

double est_period_beta(double pre_beta) {return 0.3 + 0.7 * pre_beta;}  //Relation between pre and post beta


void FindPortfolios(string incr, string folderName, string dataFolderName, Matrix Rs)
{
    //todo: Creating the BAB portfolio as if it was a stock (ETF) with daily returns

    //Define common data paths
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    //Load data from proccessed files
    Intrix iPeriod = Load_Intrix(Proccessed_FilePath_incr + "iPeriods.txt" ,-1);
    Intrix DR_iDates = Load_Intrix(Proccessed_FilePath_incr + "Rs_iDates.txt" ,-1);
    Intor Dates = Load_Intor(Proccessed_FilePath_incr + "DateList.txt");

    //Load data from prePost method
    vector <string> infoFiles = {"/year.txt", "/beta.txt", "/PERMNO.txt"};
    Tensor3 Data = Load_prePost(incr+"/"+dataFolderName, infoFiles);

    //Initiate readability variables
    int year = 0, beta = 1, PERMNO = 2;
    int pre = 0, post = 0;
    int low = 0, high = 1, all=2;
    double betaLow = 0.6, betaHigh = 0.8;

    //Initiate data variables
    int yrMin = Data[post][year][0];
    int yrMax = Data[post][year][Data[post][year].size()-1];
    vector<Intrix> lowHighID(yrMax-yrMin+1, Intrix(3, Intor(0)));
    Matrix betaLowHigh(yrMax-yrMin+1, Vector(2,0));
    Vector (yrMax-yrMin+1);


    //Calculate which stocks should go in each (low/high beta) portfolio
    int st = 0;
    for (int yr = yrMin; yr <= yrMax; ++yr) //<=max , todo: test this
    {
        int iyr = yr-yrMin;
        while(Data[post][year][st] == yr)
        {
            //Add stock to portfolio with all stocks regardless of beta
            lowHighID[iyr][all].push_back(Data[pre][PERMNO][st]);

            //calculate beta
            double hist_Beta = Data[pre][beta][st];

            //Low beta portefolio
            if(hist_Beta < betaLow)
            {
                lowHighID[iyr][low].push_back(Data[pre][PERMNO][st]);
                betaLowHigh[iyr][low] += est_period_beta(hist_Beta);
            }
            //High beta portefolio
            else if(hist_Beta > betaHigh)
            {
                lowHighID[iyr][high].push_back(Data[pre][PERMNO][st]);
                betaLowHigh[iyr][high] += est_period_beta(hist_Beta);
            }
            st++;   //next stock
        }
        //Average beta in each portefolio
        betaLowHigh[iyr][low] /= lowHighID[iyr][low].size();
        betaLowHigh[iyr][high] /= lowHighID[iyr][high].size();
    }
    //Initialising data structures
    Vector lowPortfolio(0), highPortfolio(0), allPortfolio(0), lowPortBeta(0), highPortBeta(0);
    Intor dates(0);

    //Creating portfolios
    for (int iyr = 1; iyr <= yrMax-yrMin+1; ++iyr)
    {
        //Calculating each portfolios Dly/Mly returns, consisting of earlier picked stocks
        Vector lowPort  = PortfolioReturns(lowHighID[iyr-1][low ], DR_iDates, iPeriod[iyr], Rs);
        Vector highPort = PortfolioReturns(lowHighID[iyr-1][high], DR_iDates, iPeriod[iyr], Rs);
        Vector allPort  = PortfolioReturns(lowHighID[iyr-1][all ], DR_iDates, iPeriod[iyr], Rs);

        //Saving the year of returns to the full 1927-now portfolio
        push_back(lowPortfolio, lowPort);
        push_back(highPortfolio, highPort);
        push_back(allPortfolio, allPort);

        //Calculating period length
        int periodLength = iPeriod[iyr][1]-iPeriod[iyr][0]+1;

        //Creating and saving avg pre beta in low/high portfolio
        Vector v(periodLength, betaLowHigh[iyr-1][low]);
        push_back(lowPortBeta, v);
        v = Vector(periodLength, betaLowHigh[iyr-1][high]);
        push_back(highPortBeta, v);

        //Creating and saving the date for each return
        Intor yr_dates = dateRange(iPeriod[iyr], Dates);
        push_back(dates, yr_dates);
    }

    //Create Dirs
    mkdir("Data/Output");
    mkdir("Data/Output/BAB");
    mkdir(("Data/Output/BAB/"+incr).c_str());
    string save_FilePath = "Data/Output/BAB/"+incr+"/"+folderName;
    std::filesystem::remove_all(save_FilePath.c_str());
    mkdir(save_FilePath.c_str());

    //Combine to one dataset and save
    Matrix DataSet_toSave = {highPortfolio, lowPortfolio, allPortfolio, lowPortBeta, highPortBeta, int_to_double(dates)};
    vector<string> headers = {"/highReturn.txt", "/lowReturn.txt", "/allReturn.txt", "/lowPortBeta.txt", "/highPortBeta.txt","/dates.txt"};
    SaveToCSV_transposed(save_FilePath+"/BAB_CSV.txt", DataSet_toSave, headers);
    //for (int i = 0; i<headers.size(); ++i)  Save(save_FilePath + headers[i], DataSet_toSave[i]); //Save individual files
    return;
}