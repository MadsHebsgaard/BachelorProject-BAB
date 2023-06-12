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

Matrix Calculate_Performance_for_beta(Matrix Rs, Matrix MC, Intrix iDates, const Vector& sp500, const Vector& riskFree,
        const Intor& Dates, const Intor& iPeriod, double i_Inflation, int minTradingDays, bool prePeriodNoBiasRisk, bool logarithm)
{
    Matrix Data(10, Vector(0));
    int iEnd_sp500, iStart_Rs, iLength, Active_days;
    double i_Beta, i_Alpha, i_sp500_r, i_return, i_riskFree_r, i_MrkCap;
    Vector akk_i_returns;
    int emptyCount;
    int Active_TradingDaysBeforePeriod;

    int yr = Dates[iPeriod[0]]/10000;

    for (int i = 0; i < Rs.size(); ++i)
    {
        //Krav for at aktien skal være med i periodens udregning (dvs. købes)
        if(iDates[i][1] > iPeriod[0] - minTradingDays)    continue;

        //Aktiens og sp500's start, slut og løbetid
        iStart_Rs = 1 + (iPeriod[0] - iDates[i][1]);       //index where period start in Rs[i]
        iEnd_sp500 = min(iDates[i][2], iPeriod[1]);        //sp500 end index
        iLength = iEnd_sp500 - iPeriod[0] + 1;           //iDates total

        //Dage og år aktien er aktiv (har returns != -2)
        emptyCount=0;
        for(int j = iStart_Rs; j < iStart_Rs+iLength; ++j)
            if(Rs[i][j] == -2) //Add == 0 maybe
                emptyCount++;
        Active_days = iLength - emptyCount;
        if(Active_days < 1) continue;

        //1 = all trading days required, 0=no trading days required (only when prePeriodNoBiasRisk=true, thus when backtesting pre period)
        double AccptRatio = 2.0/3.0;
        if(prePeriodNoBiasRisk)
            if(Active_days < (double) iLength * AccptRatio)  continue;

        //Active_TradingDaysBeforePeriod
        Active_TradingDaysBeforePeriod=0;
        for(int j = iStart_Rs-minTradingDays; j < iStart_Rs; ++j)
            if(Rs[i][j] > -1.5 && Rs[i][j] != 0.0)
                Active_TradingDaysBeforePeriod++;

        //Yderligere krav om at aktien er blevet handlet før tid. Hvis aktien opfylder if-statementet bliver den sprunget over.
        //Tidligere var det ATDbp <= (mTD/3.0)-1 Nu er kravet størrer
        double Acceptable_preTD_ratio = 2.0/3;  //1 = all trading days before period required, 0=no trading days B.P. required
        if(Active_TradingDaysBeforePeriod < minTradingDays * Acceptable_preTD_ratio)   continue;

        //Find Market Cap
        int yr_between = yr - round(MC[i][1]);
        size_t n_MC = MC[i].size();

        if(n_MC > yr_between+2)
            i_MrkCap = MC[i][yr_between+2];
        else if(n_MC > 2)
            i_MrkCap = MC[i][n_MC-1];
        else
            continue;

        //Stock is now safe for period and data is ready to be collected
        //Calculate akk. returns
        akk_i_returns = Calculate_akk_r(Rs[i], sp500, riskFree, iLength, iStart_Rs, iPeriod[0]);
        i_return = akk_i_returns[0];
        i_sp500_r = akk_i_returns[1];
        i_riskFree_r = akk_i_returns[2];

        //Calculate beta    //Starts minTradingDays before period
        i_Beta = Calculate_Beta(Rs[i], sp500, iPeriod[0]-minTradingDays, iEnd_sp500, iStart_Rs-minTradingDays, logarithm);

        //Calculate alpha
        i_Alpha = Calculate_Alpha(i_Beta, i_return, i_sp500_r, i_riskFree_r);

        //Save stock info
        Vector i_inf = {i_Beta, i_Alpha, i_return, Rs[i][0], i_sp500_r, i_riskFree_r, i_MrkCap, i_Inflation, (double) yr, (double) Active_days};
        for (int j = 0; j<Data.size(); ++j)
            Data[j].push_back(i_inf[j]);
    }

    return Data;
}
vector<Matrix> PrePost_Performance_for_beta(const Matrix& DR, const Matrix& MC, const Intrix& iDates, const Vector& sp500,
        const Vector& riskFree, const Intor& Dates, Intrix iPeriod, const Vector& Inflation_factor, bool logarithm)
{
    vector<Matrix> A(0);
    //size of BetaOverlap is a tradeoff of bias vs robust, makes Post_beta more robust at the cost of making Post_beta biased towards Pre_beta. Tradeoff is probably not bad.
    int BetaPrePeriod = 500, BetaOverlap = 0;  //Dly
    if (iPeriod[1][0] - iPeriod[0][0] == 12)    //Mly
    {
        BetaPrePeriod = 24;
        BetaOverlap = 0;
    }

    Matrix Pre_Data = Calculate_Performance_for_beta(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[0],
            Inflation_factor[0], BetaPrePeriod, true, logarithm);
    Matrix Post_Data = Calculate_Performance_for_beta(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[1], Inflation_factor[1],
            BetaOverlap, false, logarithm);
    return Overlapping_ID_Matrix(Pre_Data, Post_Data, 3);
}
vector<Matrix> PrePost_Calculations_for_beta(const Matrix& Rs, const Matrix& MC, const Intrix& iDates, const Vector& sp500,
        const Vector& riskFree, const Intor& Dates, const Intrix& iPeriods, const Vector& Inflation_Factor, bool logarithm)
{
    //Calculate twoPeriod_Data
    vector<Matrix> fullData(vector<Matrix>(2,Matrix(0,Vector(0))));
    for (int Period = 1; Period < iPeriods.size(); Period++)
    {
        //Calculating pre and post data
        vector<Matrix> Yly_Data = PrePost_Performance_for_beta(Rs, MC, iDates, sp500, riskFree, Dates,
                {iPeriods[Period-1], iPeriods[Period]}, {Inflation_Factor[Period-1], Inflation_Factor[Period]}, logarithm);
        push_back(fullData, Yly_Data);
    }
    return fullData;
}
void PrePost_run_for_beta(const string& incr, string folderName, Matrix Rs, bool logarithm)
{
    //Loading the necessary data
    auto start = std::chrono::system_clock::now();
    string methodName = "PrePost_for_beta";
    vector<string> logMessage, fileNames;
    Matrix MC;
    Vector Inflation_Factor, sp500, riskFree;
    Intrix iPeriods, iDates;
    Intor Dates;
    Load_Data(Rs, logMessage, MC, Inflation_Factor, iDates, sp500, riskFree, Dates, iPeriods, folderName, methodName, fileNames, incr);
    fileNames.push_back("Active_TradingDays_in_period");

    //Calculate output data in twoPeriod_Data
    vector<Matrix> fullData = PrePost_Calculations_for_beta(Rs, MC, iDates, sp500, riskFree, Dates, iPeriods, Inflation_Factor, logarithm);


    //Saving CSV and logfile
    vector<string> prePost_str = {"_pre", "_post"};
    Save_TwoDataSet_CSV(folderName+"/"+methodName+"_CSV.txt", fullData, fileNames, prePost_str);
    save_logfile(start, logMessage, folderName, logarithm);
}
