#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()
#include "load_output.h"
#include <chrono>
#include <thread>

#include "ui.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;

#pragma once

//0.1 + 0.9 * pre_beta
//0.3 + 0.7 * pre_beta
//double est_period_beta(double pre_beta) {return 0.1 + 0.9 * pre_beta;}  //Relation between pre and post beta


pair<vector<Intrix>, Tensor3> StockSelector(vector<Matrix> Data, Vector betaCond)
{
    //Initiate readability variables
    size_t year = 0, beta = 1, PERMNO = 2;
    size_t pre = 0, post = 1;
    size_t low = 0, high = 1, all=2;

    //Initiate data variables
    size_t n_LHA = 3;
    size_t yrMin = round(Data[post][year].front());
    size_t yrMax = round(Data[post][year].back());
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Tensor3 beta_LHA(yrMax-yrMin+1, Matrix(n_LHA, Vector(0)));    //Low, High and All beta portfolio

    size_t st = 0;
    for (size_t yr = yrMin; yr <= yrMax; ++yr)
    {
        size_t iyr = yr-yrMin;

        //Find middle condition for inclusion
        double beta_low_max, beta_high_min;
        if(betaCond.size()<3)
        {
            //Find median beta from year
            int s = 0;
            Vector yly_beta(0);
            while(round(Data[post][year][st+s]) == yr)
            {
                yly_beta.push_back(Data[pre][beta][st+s]);
                s++;
            }
            double beta_median = calculateMedianInRange(yly_beta, betaCond.front(), betaCond.back());
            //cout << beta_median << endl;
            beta_low_max  = beta_median;
            beta_high_min = beta_median;
        }
        else
        {
            beta_low_max  = betaCond[1];
            if(betaCond.size() == 3)    beta_high_min = betaCond[1];
            else                        beta_high_min = betaCond[2];
        }
        while(round(Data[post][year][st]) == yr)
        {
            //Finding stock id and beta
            double hist_Beta = Data[pre][beta][st];
            int stockID = round(Data[pre][PERMNO][st]);

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(stockID);
            beta_LHA[iyr][all].push_back(hist_Beta);

            //Low beta portefolio
            if(hist_Beta > betaCond.front() && hist_Beta < beta_low_max){
                ID_LHA[iyr][low].push_back(stockID);
                beta_LHA[iyr][low].push_back(hist_Beta);
            }
            //High beta portefolio
            else if(hist_Beta > beta_high_min && hist_Beta < betaCond.back()){
                ID_LHA[iyr][high].push_back(stockID);
                beta_LHA[iyr][high].push_back(hist_Beta);
            }
            st++;   //next stock
        }
    }
    return {ID_LHA, beta_LHA};
}

pair<vector<Intrix>, Tensor3> StockSelector_seeking(vector<Matrix> Data, vector<Matrix> Data_seeking, Vector info_Cond, string effect)
{
    //Initiate readability variables
    size_t year = 0, beta = 1, PERMNO = 2;
    size_t mrkCap=0, infl=1, stock_akk=2, sp500_akk=3;   //seeking data
    size_t pre = 0, post = 1;
    size_t low = 0, high = 1, all=2;

    //Initiate data variables
    size_t n_LHA = 3;
    size_t yrMin = round(Data[post][year].front());
    size_t yrMax = round(Data[post][year].back());
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Tensor3 info_LHA(yrMax-yrMin+1, Matrix(n_LHA, Vector(0)));    //Low, High and All beta portfolio

    size_t st = 0;
    for (size_t yr = yrMin; yr <= yrMax; ++yr)
    {
        size_t iyr = yr-yrMin;

        //Find middle condition for inclusion
        double low_max, high_min;

        //Find median beta from year
        int s = 0;
        Vector yly_info(0);
        while(round(Data[post][year][st+s]) == yr)
        {
            if(effect == "MC")              yly_info.push_back(log10(Data_seeking[post][mrkCap][st+s]));
            else if(effect == "momentum")   yly_info.push_back(Data_seeking[pre][stock_akk][st+s]-Data_seeking[pre][sp500_akk][st+s]);
            else cout << "Choose either 'MC' or 'momentum'";
            s++;
        }
        double info_break = calculateQuantileInRange(yly_info, info_Cond.front(), info_Cond.back(), info_Cond[1]);
        low_max  = info_break;
        high_min = info_break;
        //cout << info_break << endl;

        while(round(Data[post][year][st]) == yr)
        {
            //calculate hist_info'
            int stockID = round(Data[pre][PERMNO][st]);
            double hist_info;
            if(effect == "MC")              hist_info = log10(Data_seeking[post][mrkCap][st]);
            else if(effect == "momentum")   hist_info = Data_seeking[pre][stock_akk][st]-Data_seeking[pre][sp500_akk][st];
            else {hist_info = 0; cout << "neither 'MC' nor 'momentum' was choosen. Choose one of those.";}

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(stockID);
            info_LHA[iyr][all].push_back(hist_info);

            //Low beta portefolio
            if(info_Cond.front() < hist_info && hist_info < low_max)
            {
                ID_LHA[iyr][low].push_back(stockID);
                info_LHA[iyr][low].push_back(hist_info);
            }
                //High beta portefolio
            else if(high_min < hist_info && hist_info < info_Cond.back())
            {
                ID_LHA[iyr][high].push_back(stockID);
                info_LHA[iyr][high].push_back(hist_info);
            }
            //else    cout << "\nhist_info = " << hist_info;
            st++;   //next stock
        }
    }
    return {ID_LHA, info_LHA};
}

pair<vector<Intrix>, vector<Matrix>> StockSelector_Advanced(vector<Matrix> Data, vector<Matrix> Data_seeking, Vector betaCond)
{
    //Initiate readability variables
    size_t year = 0, beta = 1, PERMNO = 2;
    size_t mrkCap=0, infl=1, moment=2, sp500_akk=3;   //seeking data
    size_t pre = 0, post = 1;
    size_t low = 0, high = 1, all=2;

    //Initiate data variables
    size_t n_LHA = 3;
    size_t yrMin = round(Data[post][year].front());
    size_t yrMax = round(Data[post][year].back());
    vector<Intrix> ID_LHA(yrMax-yrMin+1, Intrix(n_LHA, Intor(0)));
    Tensor3 beta_LHA(yrMax-yrMin+1, Matrix(n_LHA, Vector(0)));    //Low, High and All beta portfolio

    size_t st = 0;
    for (size_t yr = yrMin; yr <= yrMax; ++yr)
    {
        size_t iyr = yr-yrMin;

        //Find middle condition for inclusion
        double beta_low_max, beta_high_min;
        if(betaCond.size()<3)
        {
            //Find median beta from year
            int s = 0;
            Vector yly_beta(0);
            while(round(Data[post][year][st+s]) == yr)
            {
                yly_beta.push_back(Data[pre][beta][st+s]);
                s++;
            }
            double beta_median = calculateMedianInRange(yly_beta, betaCond.front(), betaCond.back());
            //cout << beta_median << endl;
            beta_low_max  = beta_median;
            beta_high_min = beta_median;
        }
        else
        {
            beta_low_max  = betaCond[1];
            if(betaCond.size() == 3)    beta_high_min = betaCond[1];
            else                        beta_high_min = betaCond[2];
        }

        while(round(Data[post][year][st]) == yr)
        {
            //calculate beta
            double hist_Beta = Data[pre][beta][st];
            int stockID = round(Data[pre][PERMNO][st]);

            //Score
            double momentum = Data_seeking[pre][moment][st];
            double sp500_accumulated = Data_seeking[pre][sp500_akk][st];
            double marketCap = Data_seeking[post][mrkCap][st];
            double inflationFactor = Data_seeking[post][infl][st];
            double score = (1.5+(momentum-sp500_accumulated)/2.0)*(3.0-log10(marketCap*inflationFactor)/7.0)*(2.0-hist_Beta/6.0);

            //Add stock to portfolio with all stocks regardless of beta
            ID_LHA[iyr][all].push_back(stockID);
            beta_LHA[iyr][all].push_back(hist_Beta);

            //Low beta portefolio
            if(hist_Beta > betaCond.front() && hist_Beta < beta_low_max)
            {
                if(5.8 <= score)    //75k
                //if(6.3 <= score)
                    {
                    ID_LHA[iyr][low].push_back(stockID);
                    beta_LHA[iyr][low].push_back(hist_Beta);
                }
            }
                //High beta portefolio
            else if(hist_Beta > beta_high_min && hist_Beta < betaCond.back())
            {
                //if(3.7 < score && score < 5.2) {  //Bedst
                if(score < 5.4) {   //110k
                    ID_LHA[iyr][high].push_back(stockID);
                    beta_LHA[iyr][high].push_back(hist_Beta);
                }
            }
            st++;   //next stock
        }
    }
    return {ID_LHA, beta_LHA};
}
Matrix PortfolioCreator(const pair <vector<Intrix>, const Tensor3>& portfolio_inf, const Matrix& Rs,
        const Vector& riskFree, Intrix iPeriods, const Intrix& iDates, const Intor& Dates,
        size_t post_yrMin, size_t post_yrMax, const string& method, const bool beta_weighted)
{
    vector<Intrix> ID_LHA = portfolio_inf.first;
    Tensor3 beta_LHA = portfolio_inf.second;

    //Creating portfolios //add new function
    size_t n_LHA = ID_LHA[0].size();
    size_t iyr_start = post_yrMin-1926;
    Matrix LHA_PortData(2*n_LHA+1, Vector(0));

    for (size_t iyr = iyr_start; iyr <= post_yrMax-1926; ++iyr)  //use yrTotal instead
    {
        //Calculating period length
        size_t periodLength = iPeriods[iyr][1]-iPeriods[iyr][0]+1;

        //Calculating weigthes for low and high
        Matrix weightes(0);
        for (int LHA = 0; LHA<2; ++LHA)
        {
            //if(beta_weighted)  weightes.push_back(weightes_from_extreme(beta_LHA[iyr-iyr_start][LHA],1));
            if(beta_weighted)   weightes.push_back(RankByExtreme(beta_LHA[iyr-iyr_start][LHA], 15-30*LHA));
            else                weightes.push_back(Vector(beta_LHA[iyr-iyr_start][LHA].size(), 1));
        }

        //Equal weights for all portfolio
        weightes.push_back(Vector(beta_LHA[iyr-iyr_start][2].size(),1));

        //Calculating returns
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA], PortfolioReturns_method(
                    ID_LHA[iyr-iyr_start][LHA], weightes[LHA], iDates, iPeriods[iyr], Rs, riskFree, method));

        //Calculating betas
        for (int LHA = 0; LHA<n_LHA; ++LHA)
            push_back(LHA_PortData[LHA+n_LHA], Vector(
                    periodLength, weighted_average(beta_LHA[iyr-iyr_start][LHA], weightes[LHA])));

        //Calculating dates
        push_back(LHA_PortData[2*n_LHA], int_to_double(dateRange(iPeriods[iyr], Dates)));

        //Calculating beta_pre variance
        //double variance = Var_log(beta_LHA[iyr-iyr_start][2], 0, beta_LHA[iyr-iyr_start][2].size()-1);
        //push_back(LHA_PortData[2*n_LHA+1], Vector(periodLength, variance));
    }
    return LHA_PortData;
}
void Add_Aditinals(Matrix& LHA_PortData, vector<string>& headers, int iStart, const Matrix& Add_Vectors, const vector<string>& Add_headers)
{
    //Add additional vectors from period start
    for (auto& v: Add_Vectors)
        LHA_PortData.emplace_back(skip_First_X(v, iStart));

    //Add additional headers
    for (auto& str: Add_headers)
        headers.emplace_back(str);
}


void BackTest_run(const string& incr, string folderName, const string& dataFolderName, Matrix Rs,
        const Vector& betaCond, const string& method, const bool& beta_weighted)
{
    //std::this_thread::sleep_for(std::chrono::seconds(5)); //If needed for synchronising threads, (easy workaround)
    auto start = std::chrono::system_clock::now();

    //Loading data
    string methodName = "BackTest";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);
    logMessage.push_back("method = " + method);
    logMessage.push_back("betaCond = " + vectorToString(betaCond));
    logMessage.push_back("beta_weighted = " + to_string(beta_weighted));

    //Load data from prePost method
    vector <string> infoFiles = {"/year.txt", "/beta.txt", "/PERMNO.txt"};
    Tensor3 Data = Load_prePost(incr+"/"+dataFolderName, infoFiles);
    size_t post_yrMax = round(Data[1][0].back());
    size_t post_yrMin = round(Data[1][0].front());

    //Selecting and creating portfolios
    pair portfolio_inf = StockSelector(Data, betaCond);
    Matrix LHA_PortData = PortfolioCreator(portfolio_inf, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin,
            post_yrMax, method, beta_weighted);

    //Headers for all the information
    vector<string> headers = {"lowReturn", "highReturn", "allReturn", "lowPortBeta_pre", "highPortBeta_pre", "allPortBeta_pre", "dates"};

    //Add additional pre-known timeseries (s500, riskFree)
    Add_Aditinals(LHA_PortData, headers, Find_iDate(round(LHA_PortData[6][0]), Dates), {sp500, riskFree}, {"sp500", "riskFree"});

    //Saving CSV
    Save_CSV(folderName+"/BackTest" + "_" + method + "_CSV.txt", LHA_PortData, headers);
    bool logarithm = checkLineInFile("Data/Output/PrePost/"+incr+"/"+dataFolderName+"saveLogFile.txt", "Using log: 1");
    save_logfile(start, logMessage, folderName, logarithm);
}

void BackTest_exotic(const string& incr, string folderName, const string& dataFolderName, Matrix Rs,
        const Vector& Cond, const string& method, string effect, bool weighted)
{
    //std::this_thread::sleep_for(std::chrono::seconds(5)); //If needed for synchronising threads, (easy workaround)
    auto start = std::chrono::system_clock::now();

    //Loading data
    string methodName = "BackTest";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);
    logMessage.push_back("method = " + method);
    logMessage.push_back("Cond = " + vectorToString(Cond));

    //Load data from prePost method
    vector <string> infoFiles = {"/year.txt", "/beta.txt", "/PERMNO.txt"};
    Tensor3 Data = Load_prePost(incr+"/"+dataFolderName, infoFiles);

    size_t post_yrMax = round(Data[1][0].back());
    size_t post_yrMin = round(Data[1][0].front());

    //Advanced version utalizing BAB, small cap and momentum effects
    vector <string> infoFiles_Advanced = {"/MarketCap.txt", "/infl_factor.txt", "/akk_return.txt", "/akk_sp500.txt"};
    Tensor3 Data_Advanced = Load_prePost(incr+"/"+dataFolderName, infoFiles_Advanced);

    pair<vector<Intrix>, Tensor3> portfolio_inf;

    if (effect == "score")  portfolio_inf = StockSelector_Advanced(Data, Data_Advanced, Cond);   //MC strategy
    else                    portfolio_inf = StockSelector_seeking(Data, Data_Advanced, Cond, effect);   //score (multible)

    //Create portfolios
    Matrix LHA_PortData = PortfolioCreator(portfolio_inf, Rs, riskFree, iPeriods, iDates, Dates, post_yrMin,
            post_yrMax, method, weighted);

    //Headers for all the information
    vector<string> headers = {"lowReturn", "highReturn", "allReturn", "lowPortBeta_pre", "highPortBeta_pre", "allPortBeta_pre", "dates"};

    //Add additional pre-known timeseries (s500, riskFree)
    Add_Aditinals(LHA_PortData, headers, Find_iDate(round(LHA_PortData[6][0]), Dates), {sp500, riskFree}, {"sp500", "riskFree"});

    //Saving CSV
    Save_CSV(folderName+"/BackTest" + "_" + method + "_CSV.txt", LHA_PortData, headers);
    bool logarithm = checkLineInFile("Data/Output/PrePost/"+incr+"/"+dataFolderName+"/saveLogFile.txt", "Using log: 1");
    save_logfile(start, logMessage, folderName, logarithm);
}