#pragma once

#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "load.h"
#include <chrono>
#include <ctime>
#include <functional>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

double sum(const Vector& v);
vector<Matrix> Overlapping_ID_Matrix(Matrix A, Matrix B, int ID_row);


//dir for data
void defineFilePaths(string& incr, string& Exo_FilePath, string& Proccessed_FilePath, string& Proccessed_FilePath_inc, string& Proccessed_Dly, string& Proccessed_Mly, string& Proccessed_Yly)
{
    Exo_FilePath = "Data/Input/Exo_Files/";
    Proccessed_FilePath = "Data/Input/Processed_Files/";
    //Proccessed_FilePath = "Data/Input/Processed_Files_new2/";
    //Proccessed_FilePath = "Data/Input/Processed_Files_new3/";

    Proccessed_Dly = Proccessed_FilePath+"Dly/";
    Proccessed_Mly = Proccessed_FilePath+"Mly/";
    Proccessed_Yly = Proccessed_FilePath+"Yly/";
    Proccessed_FilePath_inc = Proccessed_FilePath + incr + "/";
}


//Core calculations
double CovSP500(Vector Stock, Vector sp500, int iStart_sp500, int iEnd_sp500, int iStart_DR)
{
    int d_iDates = iEnd_sp500 - iStart_sp500 + 1;
    double Mean_Stock = 0, Mean_sp500 = 0, Covv = 0;
    int trading_Days = 0;

    //Mean
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i+iStart_DR] != -2)
        {
            trading_Days++;
            Mean_Stock += Stock[i + iStart_DR];
            Mean_sp500 += sp500[i + iStart_sp500];
        }
    }
    Mean_Stock /= trading_Days;
    Mean_sp500 /= trading_Days;

    //CovSP500
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i+iStart_DR] != -2)
            Covv += (Stock[i+iStart_DR] - Mean_Stock) * (sp500[i + iStart_sp500] - Mean_sp500);
    }
    return Covv /(trading_Days - 1);
}
double Var(Vector sp500, int Start_number, int End_number)
{
    double Mean = 0, Variance = 0;
    //int d_dates = End_number-Start_number+1;
    //if(sp500.size() < d_dates) d_dates = sp500.size();

    int trading_Days=0;

    //Mean
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
        {
            trading_Days++;
            Mean += sp500[i];
        }
    }
    Mean /= trading_Days;

    //Variance
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
            Variance += (sp500[i] - Mean) * (sp500[i] - Mean);
    }
    return Variance /(trading_Days - 1);
}
Vector Calculate_akk_r(Vector stock, Vector sp500, Vector riskFree, int iLength, int iStart_DR, int iStart_sp500)
{
    double stockReturn=1, sp500Return=1, riskFreeReturn=1;
    int Active_days = 0;
    for (int i = 0; i < iLength; ++i)
    {
        if(stock[i+iStart_DR] < -1.5) continue;
        Active_days++;
        stockReturn *= (1 + stock[i+iStart_DR]);
        sp500Return *= (1 + sp500[i+iStart_sp500]);
        riskFreeReturn *= (1 + riskFree[i+iStart_sp500]);
    }
    return {stockReturn-1,sp500Return-1,riskFreeReturn-1};
    //return pow(stockReturn,263.5104/Active_days)-1;
}
double Calculate_Beta(const Vector& stock, const Vector& sp500, int iStart, int iEnd, int iStart_Rs)
{
    double cov = CovSP500(std::move(stock), sp500, iStart, iEnd, iStart_Rs);
    double var = Var(sp500, iStart, iEnd);
    return (cov/var);
}
double Calculate_Alpha(double beta, double mu_stock, double mu_sp500, double mu_rf) {
    return mu_stock - (mu_rf + beta * (mu_sp500 - mu_rf));
}
Matrix Calculate_Performance(Matrix Rs, Matrix MC, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, const Intor& iPeriod, double i_Inflation, int minTradingDays, bool prePeriodNoBiasRisk)
{
    Matrix Data(9, Vector(0));
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

        double AccptRatio = 2.0/3.0;  //1 = all trading days required, 0=no trading days required (only when prePeriodNoBiasRisk=true, thus when backtesting pre period)
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

        //Stock is now safe and data is collected

        //Calculate akk. returns
        akk_i_returns = Calculate_akk_r(Rs[i], sp500, riskFree, iLength, iStart_Rs, iPeriod[0]);
        i_return = akk_i_returns[0];
        i_sp500_r = akk_i_returns[1];
        i_riskFree_r = akk_i_returns[2];

        //Calculate beta    //Starts minTradingDays before period
        i_Beta = Calculate_Beta(Rs[i], sp500, iPeriod[0]-minTradingDays, iEnd_sp500, iStart_Rs-minTradingDays);

        //Calculate alpha
        i_Alpha = Calculate_Alpha(i_Beta, i_return, i_sp500_r, i_riskFree_r);

        //Save stock info
        Vector i_inf = {i_Beta, i_Alpha, i_return, Rs[i][0], i_sp500_r, i_riskFree_r, i_MrkCap, i_Inflation, (double) yr};
        for (int j = 0; j<Data.size(); ++j)
            Data[j].push_back(i_inf[j]);
    }

    return Data;
}
vector<Matrix> PrePost_Performance(const Matrix& DR, const Matrix& MC, const Intrix& iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, Intrix iPeriod, const Vector& Inflation_factor)
{
    vector<Matrix> A(0);
    //size of BetaOverlap is a tradeoff of bias vs robust, makes Post_beta more robust at the cost of making Post_beta biased towards Pre_beta. Tradeoff is probably not bad.
    int BetaPrePeriod = 500, BetaOverlap = 250;  //Dly
    if (iPeriod[1][0] - iPeriod[0][0] == 12)    //Mly
    {
        BetaPrePeriod = 24;
        BetaOverlap = 12;
    }
    //BetaPrePeriod = 0;
    //BetaOverlap = 10;


    Matrix Pre_Data = Calculate_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[0],
            Inflation_factor[0], BetaPrePeriod, true);
    Matrix Post_Data = Calculate_Performance(DR, MC, iDates, sp500, riskFree, Dates, iPeriod[1], Inflation_factor[1],
            BetaOverlap, false);
    return Overlapping_ID_Matrix(Pre_Data, Post_Data, 3);
}
Vector PortfolioReturns_method(Intor PERMNO, Intrix iDates, Intor iPeriod, Matrix Rs, Vector riskFree, const string& method)
{
    //Initiating data structures
    size_t N = PERMNO.size();
    Intor iID(N);

    vector<size_t> iStart(N);
    vector<size_t> iRun_Period(N);


    int periodLength = iPeriod[1] - iPeriod[0] + 1;
    Vector portfolio_return(periodLength);
    Vector portfolio_sum(periodLength);

    Vector shareprice(N,1);

    //Finding the indexes of the PERMNO's in Rs, to look stock up in Rs
    int ID_count = 0;
    for (int i = 0; i < Rs.size(); ++i)
        if(Rs[i][0] == PERMNO[ID_count])
        {
            iID[ID_count] = i;
            ID_count++;
        }

    //Finding the start and amount of trading days from start with stock data
    for (int i = 0; i < N; ++i)
    {
        int ID = iID[i];
        iStart[i] = (iPeriod[0] - iDates[ID][1]) + 1;   //index where period start in Rs[i]
        iRun_Period[i] = Rs[ID].size() - iStart[i];   //index where period end in Rs[i]
    }

    //Calculating portfolio returns for period
    for (int day = 0; day<periodLength; ++day)
    {
        if(method == "rebalance" || method=="RB")
        {
            double total_return=0;
            int totalLeft=0;
            for (int stock = 0; stock<N; ++stock) {
                if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)
                {
                    total_return += Rs[iID[stock]][iStart[stock] +day];
                    totalLeft++;
                }
            }
            portfolio_return[day] = total_return/totalLeft;
            continue;
        }
        //shareprice_last is updated
        Vector shareprice_last = shareprice;

        //shareprice is calculated
        for (int stock = 0; stock<N; ++stock)
        {
            //if nothing else, return is zero for the day
            double stock_return = 0;

            //If stock is still 'alive', and return is non empty (-2) return is from stock, otherwise return is from the risk free rate if method is selected
            if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)  stock_return = Rs[iID[stock]][iStart[stock] + day];
            else if(method == "RF" || method == "riskfree")                             stock_return = riskFree[iPeriod[0]+day];

            //Shareprice is now updated
            shareprice[stock] *= (1 + stock_return);
        }

        //Daily portfolio return is calculated based on either riskFree or distributed method when stock is missing data
        if(method == "RF" || method == "riskfree")
        {
            portfolio_return[day] = sum(shareprice)/sum(shareprice_last) -1;
        }
        else if(method == "DT" || method == "distributed")
        {
            double portPrice=0, portPrice_last=0;
            for (int stock = 0; stock<N; ++stock)
            {
                if(day < iRun_Period[stock] && Rs[iID[stock]][iStart[stock] + day] > -1.5)
                {
                    portPrice += shareprice[stock];
                    portPrice_last += shareprice_last[stock];
                }
            }
            portfolio_return[day] = portPrice/portPrice_last -1;
        }
    }
    return portfolio_return;
}


//Important functions
Matrix Edit_DR(Matrix A, double max_ratio, int minTradingDays)
{
    Matrix DR_ny(0);
    int traek, zero_total;
    int max_traek = 6;
    int replace = -2;
    int continue_amount = 0;

    for (auto & stock : A)
    {
        zero_total = 0;
        traek = 0;
        if (stock.size() < minTradingDays)   {continue_amount++; continue;}
        for (int j = 0; j < stock.size(); ++j)
        {
            if(stock[j] == 0 or stock[j] == replace)    {traek++; zero_total++;}
            else                 traek = 0;

            if      (traek == max_traek)  for(int k = 1; k <= max_traek; ++k) stock[j - max_traek + k] = replace;
            else if (traek > max_traek) stock[j] = replace;
        }
        if((zero_total+0.0) / ((double) stock.size()) < max_ratio)
            DR_ny.push_back(stock);
    }

    cout << "\ncontinue_amount = " << continue_amount;
    cout << "\nDR_ny size = " << DR_ny.size();
    cout << "\nZero ratio too high = " << A.size() - DR_ny.size() - continue_amount;
    cout << "\n\n";

    return DR_ny;
}
Intrix Remove_Missing_ID(Intrix A, Vector v)
{
    Intrix B(0);
    size_t k=0;
    for (auto & stock : A)
    {
        //if(k==v.size()-1) break;
        if(v[k] == stock[0]) {k++;   B.push_back(stock);}
    }
    return B;
}
Matrix Remove_Missing_ID(Matrix A, Vector v)
{
    Matrix B(0);
    size_t k=0;
    for (auto & stock : A)
    {
        //if(k==v.size()-1) break;
        if(abs(v[k] - stock[0])<0.1) {k++;   B.push_back(stock);}
    }
    return B;
}
vector<Matrix> Overlapping_ID_Matrix(Matrix A, Matrix B, int ID_row)
{
    //Finding maching indexes
    Intor A_Keep(0), B_Keep(0);
    int b = 0;
    for (int a = 0; a < A[0].size(); ++a)
    {
        while(A[ID_row][a] > B[ID_row][b])
        {
            b++;
            if(b >= B[ID_row].size()) break;
        }
        if(A[ID_row][a] == B[ID_row][b])
        {
            A_Keep.push_back(a);
            B_Keep.push_back(b);
            b++;
        }
    }
    //Keep only maching indexes
    vector<Matrix> C(2, Matrix(A.size(), Vector(A_Keep.size())));
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A_Keep.size(); ++j)
        {
            C[0][i][j] = A[i][A_Keep[j]];
            C[1][i][j] = B[i][B_Keep[j]];
        }
    }
    return C;
}


//Helper functions
double sum(const Vector& v) {
    double sum_v=0;
    for(auto& e:v)  sum_v += e;
    return sum_v;
}
Vector Matrix_Column(Matrix A, int j)
{
    Vector col(A.size());
    for (int i = 0; i < A.size(); ++i)  col[i] = A[i][j];
    return col;
}
Vector int_to_double(Intor intVec)
{
    vector<double> doubleVec(intVec.begin(), intVec.end());
    return doubleVec;
}
Vector skip_First_X(Vector v, int start_from)
{
    v.erase( v.begin(), v.size() > start_from ?  v.begin() + start_from : v.end() );
    return v;
}
void push_back(Matrix& First, const Matrix& Second)
{
    if (First.empty()) {
        First = Second;    return;
    }

    if(First.size()!=Second.size())
    {
        cout << "Push_back: Matrices have different number of rows\n";
        return;
    }
    size_t n_first;
    for (int i = 0; i < Second.size(); ++i) {
        n_first = First[i].size();
        First[i].resize(First[i].size() + Second[i].size());
        for (int j = 0; j < Second[i].size(); ++j)
            First[i][j + n_first] = Second[i][j];
    }
}
void push_back(vector<Matrix>& First, const vector<Matrix>& Second)
{
    if (First[0].empty()) {
        First = Second;    return;
    }

    if(First[0].size()!=Second[0].size())
    {
        cout << "Push_back: Matrices have different number of rows\n";
        return;
    }
    size_t n_first;
    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < Second[k].size(); ++i) {
            n_first = First[k][i].size();
            First[k][i].resize(First[k][i].size() + Second[k][i].size());
            for (int j = 0; j < Second[k][i].size(); ++j)
                First[k][i][j + n_first] = Second[k][i][j];
        }
    }
}
void push_back(Vector& first, Vector last)
{
    if(first.empty())   first = last;
    else    first.insert(first.end(), last.begin(), last.end());
}
Intor dateRange(Intor iperiod, Intor Dates)
{
    Intor dateVec(iperiod[1] - iperiod[0]+1);
    for (int i = 0; i <= iperiod[1]-iperiod[0]; ++i)
    {
        dateVec[i] = Dates[iperiod[0]+ i];
    }
    return dateVec;
}
void remove_substring(string& str, const string& sub_str)
{
    std::string::size_type i = str.find(sub_str);

    if (i != std::string::npos)
        str.erase(i, sub_str.length());
}
void remove_substring(vector<string>& v_str, const string& sub_str)
{
    for(auto& str:v_str)
        remove_substring(str, sub_str);
}
void remove_char(string& str, char char_to_remove)
{
    str.erase(remove(str.begin(), str.end(), char_to_remove), str.end()); //remove A from string
}
void remove_char(vector<string>& v_str, char char_to_remove)
{
    for(auto& str:v_str)
        remove_char(str, char_to_remove);
}
string formatNumber(int number)
{
    string str = to_string(number);
    size_t length = str.length();
    if (length <= 3) {
        return str;
    }
    size_t numSeparators = (length - 1) / 3;
    for (int i = 1; i <= numSeparators; i++) {
        str.insert(length - (i * 3), ".");
    }
    return str;
}
string formatDate(int dateInt)
{
    string dateStr = to_string(dateInt);
    string year = dateStr.substr(0, 4);
    string month = dateStr.substr(4, 2);
    string day = dateStr.substr(6, 2);
    return year + "/" + month + "/" + day;
}
int lastNumberOverK(Vector v, int k)
{
    for(int i=v.size()-1; i>=0; i--)
    {
        if (v[i] > k)
        {
            return i+1;
        }
    }
    return 0;
}


//Process files functions
int Mly_Find_iDate(int date)
{
    int year_after = date/10000 - 1926;
    int month_after = (date%10000)/100 - 1;
    return year_after*12+month_after;
}
int Find_iDate(int date, Intor Date_list)
{
    int index_apprx = 0;
    while(Date_list[index_apprx] < date)
    {
        if(index_apprx+100 < Date_list.size())
        {
            if(Date_list[index_apprx+100] <= date)
                index_apprx += 100;
            else index_apprx++;
        }
        else index_apprx++;
    }
    return index_apprx;
}
Intrix Mly_Dates_to_iDates(Intrix Rs_Dates, int start_j)
{
    for (auto & Date : Rs_Dates)
        for (int j = start_j; j < Date.size(); ++j)
            Date[j] = Mly_Find_iDate(Date[j]);
    return Rs_Dates;
}
Intrix Dly_Dates_to_iDates(Intrix Rs_Dates, const Intor& DateList, int start_j)
{
    for (auto & Date : Rs_Dates)
        for (int j = start_j; j < Date.size(); ++j)
            Date[j] = Find_iDate(Date[j], DateList);
    return Rs_Dates;
}
Intor Dly_to_Mly_DateList(Intor Dly_DateList)
{
    Intor Mly_DateList(0);
    for(int i=0; i<Dly_DateList.size(); i++)
        if (Dly_DateList[i] - Dly_DateList[i-1] > 40)
            Mly_DateList.push_back(Dly_DateList[i-1]);
    Mly_DateList.push_back(Dly_DateList[Dly_DateList.size()-1]);
    return Mly_DateList;
}
Intrix x_Periods(const string& Mly_Yly, Intor DateList)
{
    int Cond = 4000;  //Length if iPeriods
    if(Mly_Yly=="Mly" || Mly_Yly=="mly") Cond=40;   //Month long iPeriods
    if(Mly_Yly=="Yly" || Mly_Yly=="yly") Cond=4000;  //Year long iPeriods

    int i = 0, startDate=DateList[0];
    Intrix Periods(0, Intor(2));

    while (DateList[i] < 99999999)
    {
        if((DateList[i] - DateList[i-1]) > Cond)
        {
            Periods.push_back({startDate, DateList[i - 1]});
            startDate = DateList[i];
        }
        i++;
    }
    return Periods;
}
Intrix x_iPeriods(const string& length_Mly_Yly, const string& data_Dly_Mly, const Intor& DateList)
{
    Intrix iPeriods;
    Intrix Periods = x_Periods(length_Mly_Yly, DateList);

    if(data_Dly_Mly=="Dly" || data_Dly_Mly=="dly")
        iPeriods = Dly_Dates_to_iDates(Periods, DateList, 0);
    if(data_Dly_Mly=="Mly" || data_Dly_Mly=="mly")
        iPeriods = Mly_Dates_to_iDates(Periods, 0);

    return iPeriods;
}
Intrix Create_Mly_iPeriods(Intor DateList)
{
    int start = DateList.front()/10000; //rounds down to year
    int end = DateList.back()/10000;    //rounds down to year
    Intrix iPeriods(end-start+1, Intor(2));

    for (int i = 0; i<iPeriods.size(); ++i)
    {
        iPeriods[i][0] = i*12;
        iPeriods[i][1] = 11 + iPeriods[i][0];
    }
    return iPeriods;
}
Vector period_accumulate_of_Dly(Vector Dly_RFR, Intrix iPeriods)
{
    Vector acc_periods(0);
    acc_periods.reserve(iPeriods.size());

    for(auto& period:iPeriods) {
        double acc_period = 1;
        for (int i = period[0]; i<=period[1]; i++)
            acc_period *= (1.0 + Dly_RFR[i]);

        acc_periods.push_back(acc_period-1);
    }
    return acc_periods;
}
Vector DailyYearly_to_DailyDaily_Return(Vector DailyYearly, Intrix iPeriods)
{
    int TradingDays;
    for(auto &period:iPeriods)
    {
        TradingDays = period[1]-period[0]+1;
        for (int i = period[0]; i <= period[1]; ++i)
        {
            DailyYearly[i] = pow(1.0 + DailyYearly[i],1.0/TradingDays)-1.0;
        }
    }
    return DailyYearly;
}
Matrix MarketCap_Mly_to_Yly(Matrix MCm)
{
    Matrix MCy(0);
    int FirstMonth;
    for (auto & MCi : MCm) {
        Vector Stock(2);
        Stock[0]  = MCi[0];              //ID
        FirstMonth = fmod(MCi[1],100);   // MCm[MCi][1] % 100
        int yr = (int) MCi[1] / 100;
        Stock[1] = yr;
        if (FirstMonth != 1) Stock[1]++;

        if (FirstMonth == 1) FirstMonth = 13;
        for (int j = 2+(13-FirstMonth); j < MCi.size(); j+=12) {
            Stock.push_back(MCi[j]);
        }
        if((MCi.size()-2-(13-FirstMonth))%12 == 0)   Stock.push_back(MCi[MCi.size()-1]);
        //if(Stock.size()>2)  //Ignore stock that die before Market Cap is recorded.
        MCy.push_back(Stock);
    }
    return MCy;
}
Vector Inflation_Factors_from_yrly_inf(Vector yrly_inflation)
{
    Vector infl_Factors(yrly_inflation.size()+1);
    infl_Factors[0] = 1;
    for (int yr = 1; yr < yrly_inflation.size()+1; ++yr) {
        infl_Factors[yr] = infl_Factors[yr-1] / (1+yrly_inflation[yr-1]);
    }
    return infl_Factors;
}
Matrix combine_First_with_Mly_MC(Matrix MC, Matrix fMC)
{
    Matrix combined_MC;
    size_t i=0, j=0, n_MC=MC.size(), n_fMC=fMC.size();
    int last_ID=9999;

    while(i < n_MC-1 || j < n_fMC-1) //todo: <= or <
    {
        while  (MC[i][0] <= last_ID && i<n_MC -1)   i++;
        while (fMC[j][0] <= last_ID && j<n_fMC-1)  j++;

        if    (MC[i][0] <= fMC[j][0] && MC[i][0] > last_ID)
        {
            combined_MC.push_back(MC[i]);
            last_ID = round(MC[i][0]);
        }
        else if (fMC[j][0] > last_ID)
        {
            Vector stock(3);
            for (int k = 0; k<3; ++k)   stock[k] = fMC[j][k];
            combined_MC.push_back(stock);
            last_ID = round(stock[0]);
        }
    }
    return combined_MC;
}
Matrix Missing_MC_to_Zero(Matrix pMC2, Intrix Rs_Dates)
{
    Matrix pMC_comb(0);
    int j=0;
    for(int i=0; i<Rs_Dates.size(); i++)
    {
        if(pMC2[j][0] == Rs_Dates[i][0])
        {
            pMC_comb.push_back(pMC2[j]);
            j++;
        }
        else if(Rs_Dates[i][0] < pMC2[j][0]){
            //cout << Rs_Dates[i][0] << "<" << pMC2[j][0] << endl;
            pMC_comb.push_back({Rs_Dates[i][0]*1.0, Rs_Dates[i][1]/100.0, 0});
        }
        else if(Rs_Dates[i][0] > pMC2[j][0])
        {
            //cout << Rs_Dates[i][0] << ">" << pMC2[j][0] << endl;
            i--;    //Stay at index
            j++;    //Go one index down
        }
    }
    return pMC_comb;
}
void whatToLoad(string& Dly_Mly_Both, bool& Dly, bool& Mly)
{
    if     (Dly_Mly_Both=="Dly" || Dly_Mly_Both=="dly" )    Dly=true;
    else if(Dly_Mly_Both=="Mly" || Dly_Mly_Both=="mly" )    Mly=true;
    else if(Dly_Mly_Both=="Both"|| Dly_Mly_Both=="both")   {Dly=true; Mly=true;}
}
bool areFilesExistInDirectory(const std::vector<std::string>& filenames, const std::string& directoryPath)
{
    bool allExist = true;
    for (const auto& filename : filenames) {
        const auto filepath = directoryPath + filename;
        if (!filesystem::exists(filepath)) {
            std::cerr << "File " << filepath << " does not exist\n";
            allExist = false;
        }
    }
    if(allExist)    return true;
    else            return false;
}
int numSubdirsInDir(const string& dirPath)
{
    const std::filesystem::path dir_path = dirPath;
    int num_subdirs = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.is_directory() && entry.path().filename().string()[0] != '.') {
            ++num_subdirs;
        }
    }
    return num_subdirs;
}