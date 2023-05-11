#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()
#include "load_output.h"
#include <iostream>
#include <chrono>
#include <thread>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;

#pragma once

//0.1 + 0.9 * pre_beta
//0.3 + 0.7 * pre_beta
double est_period_beta(double pre_beta) {return 0.1 + 0.9 * pre_beta;}  //Relation between pre and post beta

pair<vector<Intrix>, Matrix> StockSelector(vector<Matrix> Data, Vector betaCond)
{
    //Initiate readability variables
    int year = 0, beta = 1, PERMNO = 2;
    int pre = 0, post = 1;
    int low = 0, high = 1, all=2;

    //Initiate data variables
    int n_LHA = 3;
    int yrMin = Data[post][year][0];
    int yrMax = Data[post][year][Data[post][year].size()-1];
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Matrix beta_LHA(yrMax-yrMin+1, Vector(n_LHA,0));    //Low, High and All beta portfolio

    int st = 0;
    for (int yr = yrMin; yr <= yrMax; ++yr)
    {
        int iyr = yr-yrMin;
        while(Data[post][year][st] == yr)
        {
            //calculate beta
            double hist_Beta = Data[pre][beta][st];

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(Data[pre][PERMNO][st]);
            beta_LHA[iyr][all] += (hist_Beta);

            //Low beta portefolio
            if(hist_Beta > betaCond[0] && hist_Beta < betaCond[1])
            {
                ID_LHA[iyr][low].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][low] += (hist_Beta);
            }
                //High beta portefolio
            else if(hist_Beta > betaCond[2] && hist_Beta < betaCond[3])
            {
                ID_LHA[iyr][high].push_back(Data[pre][PERMNO][st]);
                beta_LHA[iyr][high] += (hist_Beta);
            }
            st++;   //next stock
        }
        //Average beta in each portefolio
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            beta_LHA[iyr][LHA] /= ID_LHA[iyr][LHA].size();
    }
    return {ID_LHA, beta_LHA};
}
pair<vector<Intrix>, Matrix> StockSelector_Advanced(vector<Matrix> Data, vector<Matrix> Data_seeking, Vector betaCond)
{
    //Initiate readability variables
    int year = 0, beta = 1, PERMNO = 2;
    int mrkCap=0, infl=1, moment=2, sp500_akk=3;   //seeking data
    int pre = 0, post = 1;
    int low = 0, high = 1, all=2;

    //Initiate data variables
    int n_LHA = 3;
    int yrMin = Data[post][year][0];
    int yrMax = Data[post][year][Data[post][year].size()-1];
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Matrix beta_LHA(yrMax-yrMin+1, Vector(n_LHA,0));    //Low, High and All beta portfolio

    int st = 0;
    for (int yr = yrMin; yr <= yrMax; ++yr)
    {
        int iyr = yr-yrMin;
        while(Data[post][year][st] == yr)
        {
            //calculate beta
            double hist_Beta = Data[pre][beta][st];

            //Score
            double momentum = Data_seeking[pre][moment][st];
            double sp500_accumulated = Data_seeking[pre][sp500_akk][st];
            double marketCap = Data_seeking[post][mrkCap][st];
            double inflationFactor = Data_seeking[post][infl][st];
            double score = (1.5+(momentum-sp500_accumulated)/2.0)*(3.0-log10(marketCap*inflationFactor)/7.0)*(2.0-hist_Beta/6.0);

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(Data[pre][PERMNO][st]);
            beta_LHA[iyr][all] += hist_Beta;

            //Low beta portefolio
            if(hist_Beta > betaCond[0] && hist_Beta < betaCond[1])
            {
                if(5.8 <= score)    //75k
                //if(6.3 <= score)
                    {
                    ID_LHA[iyr][low].push_back(Data[pre][PERMNO][st]);
                    beta_LHA[iyr][low] += (hist_Beta);
                }
            }
                //High beta portefolio
            else if(hist_Beta > betaCond[2] && hist_Beta < betaCond[3])
            {
                //if(3.7 < score && score < 5.2) {  //Bedst
                if(score < 5.4) {   //110k
                    ID_LHA[iyr][high].push_back(Data[pre][PERMNO][st]);
                    beta_LHA[iyr][high] += (hist_Beta);
                }
            }
            st++;   //next stock
        }
        //Average beta in each portefolio
        for (int LHA = 0; LHA<n_LHA; ++LHA)
        {
            cout << LHA << " = " << ID_LHA[iyr][LHA].size() << endl;
            beta_LHA[iyr][LHA] /= ID_LHA[iyr][LHA].size();
        }
    }
    return {ID_LHA, beta_LHA};
}
Matrix LHA_PortfolioCreator(vector<Intrix> ID_LHA, Matrix beta_LHA, Matrix Rs, Vector riskFree, Intrix iPeriods, Intrix iDates, Intor Dates, int post_yrMin, int post_yrMax, string method)
{
    //Creating portfolios //add new function
    int n_LHA = ID_LHA[0].size();
    int iyr_start = post_yrMin-1926;
    Matrix LHA_PortData(2*n_LHA+1, Vector(0));
    for (int iyr = iyr_start; iyr <= post_yrMax-1926; ++iyr)  //use yrTotal instead
    {
        //Calculating period length
        int periodLength = iPeriods[iyr][1]-iPeriods[iyr][0]+1;

        //Calculating returns
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA], PortfolioReturns_method(ID_LHA[iyr-iyr_start][LHA], iDates, iPeriods[iyr], Rs, riskFree, method));

        //Calculating betas
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA+n_LHA], Vector(periodLength, beta_LHA[iyr-iyr_start][LHA]));

        //Calculating dates
        push_back(LHA_PortData[2*n_LHA], int_to_double(dateRange(iPeriods[iyr], Dates)));
    }
    return LHA_PortData;
}
void Add_Aditinals(Matrix& LHA_PortData, vector<string>& headers, int iStart, Matrix Add_Vectors, vector<string> Add_headers)
{
    //Add additional vectors from period start
    for (auto& v: Add_Vectors)
        LHA_PortData.emplace_back(skip_First_X(v, iStart));

    //Add additional headers
    for (auto& str: Add_headers)
        headers.emplace_back(str);
}
void BackTest_run(string incr, string folderName, string dataFolderName, Matrix Rs, Vector betaCond, string method)
{
    auto start = std::chrono::system_clock::now();


    //Loading data
    string methodName = "BackTest";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);
    //logMessage.push_back(to_string(BetaCond));

    //Load data from prePost method
    vector <string> infoFiles = {"/year.txt", "/beta.txt", "/PERMNO.txt"};
    Tensor3 Data = Load_prePost(incr+"/"+dataFolderName, infoFiles);

    //Calculate which stocks should go in each (low/high beta) portfolio
    auto [ID_LHA, beta_LHA] = StockSelector(Data, betaCond);

    //Advanced version utalizing BAB, small cap and momentum effects
    //vector <string> infoFiles_Advanced = {"/MarketCap.txt", "/infl_factor.txt", "/akk_return.txt", "/akk_sp500.txt"};
    //Tensor3 Data_Advanced = Load_prePost(incr+"/"+dataFolderName, infoFiles_Advanced);
    //auto [ID_LHA, beta_LHA] = StockSelector_Advanced(Data, Data_Advanced, betaCond);   //Advanced strategi

    int post_yrMax = Data[1][0][Data[0][0].size()-1];
    int post_yrMin = Data[1][0][0];
    //int yrTotal = post_yrMax - post_yrMin + 1;

    //Calculate performance of LHA portfolios, their Beta values and the year
    Matrix LHA_PortData = LHA_PortfolioCreator(ID_LHA, beta_LHA, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin, post_yrMax, method);

    //Headers for all the information
    vector<string> headers = {"lowReturn", "highReturn", "allReturn", "lowPortBeta_pre", "highPortBeta_pre", "allPortBeta_pre", "dates"};

    //Add additional pre-known timeseries (s500, riskFree)
    Add_Aditinals(LHA_PortData, headers, Find_iDate(LHA_PortData[6][0], Dates), {sp500, riskFree}, {"sp500", "riskFree"});

    //Saving CSV
    Save_CSV(folderName+"/BackTest" + "_" + method + "_CSV.txt", LHA_PortData, headers);
    save_logfile(start, logMessage, folderName);
}