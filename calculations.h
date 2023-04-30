#pragma once

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

void push_back(vector<Matrix>& First, const vector<Matrix>& Second);
void push_back(Matrix& First, const Matrix& Second);
double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR);      //TODO
pair<double,double> Calculate_Return(Vector stock, int iLength, int iStart_DR);             //TODO
double CovSP500(Vector Stock, Vector sp500, int iStart_sp500, int iEnd_sp500, int iStart_DR);   //TODO: Fix start date and end date
double Var(Vector sp500, int Start_number, int End_number);                         //TODO: Fix start date and end date
int day_after1926(int date);
Intor Column_of_Intrix(Intrix A, int n);
double Vector_Average(const Vector& v);
Vector Matrix_Column(Matrix A, int j);
void Save(const string& fn, const Vector& v);
Intor Load_Intor(const string& fn);
Intrix Dates_to_iDates(function<int(const int&, const Intor&)> Find_iDate, Intrix Rs_Dates, Intor DateList, int start_j);
double average(Vector v)
{
    double sum = 0;
    for(auto& e:v)
        sum+=e;
    return sum/(v.size()+0.0);
}
double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR)
{
    double stockCov = CovSP500(stock, sp500, iStart, iEnd, iStart_DR);
    double stockSP500var = Var(sp500, iStart, iEnd);
    return (stockCov/stockSP500var);
}
Vector Calculate_akk_Return(Vector stock, Vector sp500, Vector riskFree, int iLength, int iStart_DR, int iStart_sp500)
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
int day_after1926(int date)
{
    int year_after = date/10000 - 1926;
    int month_after = (date%10000)/100 - 1;
    int days_after = date%100-2;
    return 263 * year_after + month_after * 23 + days_after * 13/16;
} //TODO: Check if it works decently for years after 1927
Intor Column_of_Intrix(Intrix A, int n)
{
    Intor v(A.size());
    for (int i = 0; i < A.size(); ++i) {
        v[i] = A[i][n];
    }
    return v;
}
int Dly_Find_iDate(int date, Intor Date_list)
{
    //Evt. pushback 99999999 instead of having datelist ending in 9999999
    int index_apprx = day_after1926(date);
    while(Date_list[index_apprx] < date)    index_apprx += 15;
    while(Date_list[index_apprx] > date)    index_apprx -= 1;
    return index_apprx;
}
int Mly_Find_iDate(int date)
{
    int year_after = date/10000 - 1926;
    int month_after = (date%10000)/100 - 1;
    return year_after*12+month_after;
}
double Vector_Average(const Vector& v)
{
    double sum=0;
    for (double e : v)  sum+=e;
    return sum/v.size();
}
Vector Matrix_Column(Matrix A, int j)
{
    Vector col(A.size());
    for (int i = 0; i < A.size(); ++i)  col[i] = A[i][j];
    return col;
}
Intrix Dly_Dates_to_iDates(Intrix Rs_Dates, Intor DateList, int start_j)
{
    for (auto & Date : Rs_Dates)
        for (int j = start_j; j < Date.size(); ++j)
            Date[j] = Dly_Find_iDate(Date[j], DateList);
    return Rs_Dates;
}
Intrix Mly_Dates_to_iDates(Intrix Rs_Dates, int start_j)
{
    for (auto & Date : Rs_Dates)
        for (int j = start_j; j < Date.size(); ++j)
            Date[j] = Mly_Find_iDate(Date[j]);
    return Rs_Dates;
}
/*Intrix Dly_Dates_to_iDates(int (*Find_iDate)(int, Intor) ,Intrix Rs_Dates, Intor DateList, int start_j)
{
    for (auto & Date : Rs_Dates)
    {
        for (int j = start_j; j < Date.size(); ++j)
        {
            Date[j] = Find_iDate(Date[j], DateList);
        }
    }
    return Rs_Dates;
}*/
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
        if((zero_total+0.0) / (stock.size()+0.0) < max_ratio)  {DR_ny.push_back(stock);           /*cout << endl << "zero_total:" << zero_total << "  size: " << A[stock].size();*/}
        //else cout << endl << stock << "   " << A[stock][0] << "  " << (zero_total+0.0)/A[stock].size();
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
Intrix dPeriods_From_iPeriod(Intrix iPeriod)
{
    Intor Dates = Load_Intor("Dly_DateList.txt");
    for(int i =0; i<iPeriod.size(); i++)
    {
        iPeriod[i][0] = Dates[iPeriod[i][0]];
        iPeriod[i][1] = Dates[iPeriod[i][1]];
    }
    return iPeriod;
}
Intrix x_Periods(string Mly_Yly, Intor DateList)
{
    int Cond = 4000;  //Length if iPeriods
    if(Mly_Yly=="Mly" || Mly_Yly=="mly") Cond=40;   //Month long iPeriods
    if(Mly_Yly=="Yly" || Mly_Yly=="yly") Cond=4000;  //Year long iPeriods

    int i = 0, startDate=DateList[0];
    Intrix Periods(0);

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
Intrix x_iPeriods(string length_Mly_Yly, string data_Dly_Mly, Intor DateList)
{
    Intrix iPeriods;
    Intrix Periods = x_Periods(length_Mly_Yly, DateList);

    if(data_Dly_Mly=="Dly" || data_Dly_Mly=="dly")
        iPeriods = Dly_Dates_to_iDates(Periods, DateList, 0);
    if(data_Dly_Mly=="Mly" || data_Dly_Mly=="mly")
        iPeriods = Mly_Dates_to_iDates(Periods, 0);

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
int min_Val_Index(Vector v)
{
    int minIndex = 0;
    for (int i = 1; i < v.size(); ++i)
        if(v[i] < v[minIndex])
            minIndex = i;
    return minIndex;
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
int DaysBetween(int iDate_from, int iDate_to) {
    //TODO: iDates to dDates
    int date1 = iDate_from;
    int date2 = iDate_to;

    // extract year, month, and day from date1
    int year1 = date1 / 10000;
    int month1 = (date1 / 100) % 100;
    int day1 = date1 % 100;

    // extract year, month, and day from date2
    int year2 = date2 / 10000;
    int month2 = (date2 / 100) % 100;
    int day2 = date2 % 100;

    // calculate the number of days between the two dates
    int days1 = year1 * 365 + std::floor((year1 - 1) / 4.0) - std::floor((year1 - 1) / 100.0) + std::floor((year1 - 1) / 400.0) + std::floor((367 * month1 - 362) / 12.0) + ((month1 <= 2) ? 0 : ((year1 % 4 == 0 && year1 % 100 != 0) || year1 % 400 == 0) ? -1 : -2) + day1;
    int days2 = year2 * 365 + std::floor((year2 - 1) / 4.0) - std::floor((year2 - 1) / 100.0) + std::floor((year2 - 1) / 400.0) + std::floor((367 * month2 - 362) / 12.0) + ((month2 <= 2) ? 0 : ((year2 % 4 == 0 && year2 % 100 != 0) || year2 % 400 == 0) ? -1 : -2) + day2;
    return std::abs(days1 - days2);
}
double Calculate_StockAlpha(double stockBeta, double stock_Return_akk, double sp500_Return_akk, double RiskFree_Return_akk, double Active_years)
{
    //TODO: Tilbagediskunter med risikofrie rente?
    double stock_avg_Return = pow(1 + RiskFree_Return_akk, Active_years) - 1;
    double sp500_avg_Return = pow(1 + sp500_Return_akk, 1/Active_years) - 1;
    double riskFree_avg_Return = pow(1 + RiskFree_Return_akk, 1/Active_years) - 1;

    /*//Old way
    double SML = stockBeta * (sp500_avg_Return - riskFree_avg_Return) + riskFree_avg_Return;
    if (SML > -1)
        stockAlpha = stock_Return_akk - (pow(1.0 + SML, Active_years) - 1.0);
    else    //Avoids nan. Maybe makes sence?
        stockAlpha = stock_Return_akk - (-pow(-(1.0 + SML), Active_years) - 1.0);*/

    //stockAlpha = stock_avg_Return - (stockBeta * (sp500_avg_Return - riskFree_avg_Return) + riskFree_avg_Return);

    return stock_Return_akk - (RiskFree_Return_akk + stockBeta * (sp500_Return_akk - RiskFree_Return_akk));
}
vector<Intrix> SplitPeriods(const Intrix& iPeriods, int x, bool rest_in_first_period)
{
    int numRowsPerVec = floor(iPeriods.size() / x);
    int rest = iPeriods.size() - x * numRowsPerVec;

    vector<vector<vector<int>>> threeDimVec(x, vector<vector<int>>(numRowsPerVec, vector<int>(2)));

    int currRow = 0;

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < numRowsPerVec; j++) {
            threeDimVec[i][j][0] = iPeriods[currRow][0];
            threeDimVec[i][j][1] = iPeriods[currRow][1];
            currRow++;
        }
        if(i==0)
        {
            if(rest_in_first_period)
            {
                if (rest > 0)
                {
                    for (int i = 0; i < rest; ++i)
                    {
                        threeDimVec[0].push_back({iPeriods[currRow][0],iPeriods[currRow][1]});
                        currRow++;
                    }
                }
            }
        }
    }

    cout << "Periods of years was made to be (" << x << " x " << numRowsPerVec << ")";
    if(rest > 0)
    {
        if(rest_in_first_period)
            cout << ", first period have " << numRowsPerVec+rest << " periods.\n";
        else
            cout << ", first period also have " << numRowsPerVec << " periods, so " << rest << " yeas was not included in any period.\n";
    }
    else
        cout << ".\n";

    return threeDimVec;
}
Matrix Beta_Alpha_Calculate(Matrix DR, Matrix MC, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, const Intor& iPeriod, double Inflation_factor, int minTradingDaysBeforePeriod, bool prePeriodNoBiasRisk)
{
    Vector beta(0), stock_return_akk(0), sp500_return_akk(0), alpha(0), PERMNO(0), riskFree_Return_akk(0), marketCap(0), inflation_factor(0), year_v(0);;
    int iStart_sp500, iEnd_sp500, iStart_DR, iLength, Active_days;
    double stockBeta, stockAlpha, sp500_Return_akk, stock_Return_akk, RiskFree_Return_akk, Active_years, stock_MrkCap;
    Vector akk_stock_sp500_RiskFree_Return;
    int emptyCount;
    int Active_TradingDaysBeforePeriod;

    beta.reserve(DR.size()); alpha.reserve(DR.size()); stock_return_akk.reserve(DR.size()); PERMNO.reserve(DR.size()); sp500_return_akk.reserve(DR.size()); riskFree_Return_akk.reserve(DR.size());

    for (int i = 0; i < DR.size(); ++i)
    {
        //Krav for at aktien skal være med i periodens udregning (dvs. købes)
        if(iDates[i][1] > iPeriod[0] - minTradingDaysBeforePeriod)    continue;

        //Aktiens og sp500's start, slut og løbetid
        iStart_DR = 1 + (iPeriod[0] - iDates[i][1]);       //index where period start in DR[i]
        iStart_sp500 = iPeriod[0];                         //sp500 start index      //iStart_sp500 = iDates[i][1] + iStart_DR-1;
        iEnd_sp500 = min(iDates[i][2], iPeriod[1]);        //sp500 end index
        iLength = iEnd_sp500 - iStart_sp500 + 1;           //iDates total

        //Dage og år aktien er aktiv (har returns != -2)
        emptyCount=0;
        for(int j = iStart_DR; j < iStart_DR+iLength; ++j)
            if(DR[i][j] == -2) //Add == 0 maybe
                emptyCount++;
        Active_days = iLength - emptyCount;
        if(Active_days < 1) continue;

        double Acceptable_TD_ratio = 1.0/3.0;  //1 = all trading days required, 0=no trading days required (only when prePeriodNoBiasRisk=true, thus when backtesting pre period)
        if(prePeriodNoBiasRisk)
            if(Active_days < iLength * Acceptable_TD_ratio)  continue;

        //Active_TradingDaysBeforePeriod
        Active_TradingDaysBeforePeriod=0;
        for(int j = iStart_DR-minTradingDaysBeforePeriod; j < iStart_DR; ++j)
            if(DR[i][j] > -1.5 && DR[i][j] != 0.0)
                Active_TradingDaysBeforePeriod++;

        //Yderligere krav om at aktien er blevet handlet før tid. Hvis aktien opfylder if-statementet bliver den sprunget over.
        //Tidligere var det ATDbp <= (mTD/3.0)-1 Nu er kravet størrer
        double Acceptable_preTD_ratio = 1.0/3.0;  //1 = all trading days before period required, 0=no trading days B.P. required
        if(Active_TradingDaysBeforePeriod < minTradingDaysBeforePeriod * Acceptable_preTD_ratio)   continue;


        //Calculate beta    //Starts minTradingDaysBeforePeriod before period
        stockBeta = Calculate_stockBeta(DR[i], sp500, iStart_sp500-minTradingDaysBeforePeriod, iEnd_sp500, iStart_DR-minTradingDaysBeforePeriod);

        //Calculate akk. returns
        akk_stock_sp500_RiskFree_Return = Calculate_akk_Return(DR[i], sp500, riskFree, iLength, iStart_DR, iStart_sp500);
        stock_Return_akk = akk_stock_sp500_RiskFree_Return[0];
        sp500_Return_akk = akk_stock_sp500_RiskFree_Return[1];
        RiskFree_Return_akk = akk_stock_sp500_RiskFree_Return[2];

        //Calculate alpha
        Active_years = 0; //Active_years = (DaysBetween(Dates[iStart_sp500], Dates[iStart_sp500+Active_days]) + 0.0) / 365.24;
        stockAlpha = Calculate_StockAlpha(stockBeta, stock_Return_akk, sp500_Return_akk, RiskFree_Return_akk, Active_years);

        //Find Market Cap.
        int yr = Dates[iPeriod[0]]/10000;
        int yr_between = yr - MC[i][1];
        int n_MC = MC[i].size();

        if(n_MC > yr_between+2)
            stock_MrkCap = MC[i][yr_between+2];
        else if(n_MC > 2)
            stock_MrkCap = MC[i][n_MC-1];
        else
            continue;

        if(stock_MrkCap < 100)
            cout << "Mrk Cap = " << stock_MrkCap << " For ID: " << DR[i][0] << " in year = " << yr << endl;


        //Last index if real index later then MC[i].size()
            //If last index is then <2, continue


        //Save data to vectors
        beta.push_back(stockBeta);
        alpha.push_back(stockAlpha);
        stock_return_akk.push_back(stock_Return_akk);
        sp500_return_akk.push_back(sp500_Return_akk);
        riskFree_Return_akk.push_back(RiskFree_Return_akk);
        marketCap.push_back(stock_MrkCap);
        inflation_factor.push_back(Inflation_factor);
        year_v.push_back(yr);
        PERMNO.push_back(DR[i][0]);
    }
    //double avg_sp500 = Calculate_akk_Return(sp500, sp500, riskFree, riskFree.size(), 0, 0)[1];
    //double avg_riskFree = Calculate_akk_Return(sp500, sp500, riskFree, riskFree.size(), 0, 0)[2];
    return {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk, marketCap, inflation_factor, year_v};
}
inline bool filePath_exists(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
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
bool areFilesExistInEitherDirectory(const vector<string>& filenames, const vector<string>& directoryPaths)
{
    bool allExist = true;
    for (const auto& filename : filenames) {
        for(const auto& directoryPath : directoryPaths)
        {
            const auto filepath = directoryPath + filename;
            if (!filesystem::exists(filepath)) {
                std::cerr << "File " << filepath << " does not exist\n";
                allExist = false;
            }
        }
    }
    if(allExist)    return true;
    else            return false;
}
void push_back(Matrix& First, const Matrix& Second)
{
    if (First.size()==0) {
        First = Second;    return;
    }

    if(First.size()!=Second.size())
    {
        cout << "Push_back: Matrices have different number of rows\n";
        return;
    }
    double n_first;
    for (int i = 0; i < Second.size(); ++i) {
        n_first = First[i].size();
        First[i].resize(First[i].size() + Second[i].size());
        for (int j = 0; j < Second[i].size(); ++j)
            First[i][j + n_first] = Second[i][j];
    }
    return;
}
void push_back(vector<Matrix>& First, const vector<Matrix>& Second)
{
    if (First[0].size()==0) {
        First = Second;    return;
    }

    if(First[0].size()!=Second[0].size())
    {
        cout << "Push_back: Matrices have different number of rows\n";
        return;
    }
    double n_first;
    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < Second[k].size(); ++i) {
            n_first = First[k][i].size();
            First[k][i].resize(First[k][i].size() + Second[k][i].size());
            for (int j = 0; j < Second[k][i].size(); ++j)
                First[k][i][j + n_first] = Second[k][i][j];
        }
    }
    return;
}
vector<Matrix> Overlapping_ID_Matrix_Array(Matrix A, Matrix B, int ID_row)
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
vector<Matrix> TwoPeriod_Calc(Matrix DR, Matrix MC, Intrix iDates, Vector sp500, Vector riskFree, Intor Dates, Intor Pre_iPeriod, Intor iPeriod, Vector Inflation_factor)
{
    vector<Matrix> A(0);
    //size of BetaOverlap is a tradeoff of bias vs robust, makes Post_beta more robust at the cost of making Post_beta biased towards Pre_beta. Tradeoff is probably not bad.
    int BetaOverlap = 20;                   //Dly
    if (iPeriod[0] - Pre_iPeriod[0] == 12)  //Mly
        BetaOverlap = 3;

    Matrix Pre_Period_values = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, Pre_iPeriod, Inflation_factor[0], 0, true);
    Matrix Period_values = Beta_Alpha_Calculate(DR, MC, iDates, sp500, riskFree, Dates, iPeriod, Inflation_factor[1], BetaOverlap, false);
    return Overlapping_ID_Matrix_Array(Pre_Period_values, Period_values, 3);
}
int numFilesInDir(string dirPath)
{
    const std::filesystem::path dir_path = dirPath;
    int num_files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.path().filename().string()[0] != '.') {
            ++num_files;
        }
    }
return num_files;
}
int numSubdirsInDir(string dirPath)
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
Matrix MarketCap_Mly_to_Yly(Matrix MCm)
{
    Matrix MCy(0);
    int FirstMonth;
    for (int i = 0; i < MCm.size(); ++i) {
        Vector Stock(2);
        Stock[0]  = MCm[i][0];              //ID
        FirstMonth = fmod(MCm[i][1],100);   // MCm[i][1] % 100
        int yr = MCm[i][1] / 100;
        Stock[1] = yr;
        if (FirstMonth != 1) Stock[1]++;

        if (FirstMonth == 1) FirstMonth = 13;
        for (int j = 2+(13-FirstMonth); j < MCm[i].size(); j+=12) {
            Stock.push_back(MCm[i][j]);
        }
        if((MCm[i].size()-2-(13-FirstMonth))%12 == 0)   Stock.push_back(MCm[i][MCm[i].size()-1]);
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
int countOfElements(Matrix A)
{
    int count = 0;
    for(auto v:A) {
        count += v.size();
    }
    return count;
}
string formatNumber(int number)
{
    string str = to_string(number);
    int length = str.length();
    if (length <= 3) {
        return str;
    }
    int numSeparators = (length - 1) / 3;
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
Vector PortfolioReturns(Intor PERMNO, Intrix iDates, Intor iPeriod, Matrix DR)
{
    int N = PERMNO.size();
    Vector includeID(N);
    int countID = 0;
    for (int i = 0; i < DR.size(); ++i)
    {
        if(DR[i][0] == PERMNO[countID])
        {
            includeID[countID] = i;
            countID++;
        }
    }
    Intor iStart_DR(N);
    Intor iRun_Period(N);
    //Intor iStart_Period_ID(N);          iStart_Period_ID[i] = max(iPeriod[0], iDates[ID][0]) - iPeriod[0];


    for (int i = 0; i < N; ++i) {
        int ID = includeID[i];
        iStart_DR  [i] = (iPeriod[0] - iDates[ID][1]) + 1;                //index where period start in DR[i]
        iRun_Period[i] = DR[ID].size() - iStart_DR[i];   //index where period end in DR[i]
    }
    int periodLength = iPeriod[1] - iPeriod[0] + 1;

    Vector portfolio_return(periodLength);
    Vector portfolio_kurs(periodLength);
    Vector stock_kurs(N,1);

    for (int day = 0; day<periodLength; ++day) {

        Vector stock_kurs_new(0);

        for (int stock = 0; stock<N; ++stock) {
            if(day < iRun_Period[stock])
            {
                double stock_return = DR[includeID[stock]][iStart_DR[stock] + day];
                if(stock_return > -1.5)   stock_kurs[stock] *= (1.0 + stock_return);
                stock_kurs_new.push_back(stock_kurs[stock]);
            }
        }
        portfolio_kurs[day] = average(stock_kurs_new);

        if(day != 0)    portfolio_return[day] = portfolio_kurs[day]/portfolio_kurs[day-1] -1;
        else            portfolio_return[day] = portfolio_kurs[day]/1 -1;
    }

    return portfolio_return;
}
void push_back(Vector& first, Vector& last)
{
    first.insert(first.end(), last.begin(), last.end());
}
void push_back(Intor& first, Intor& last)
{
    first.insert(first.end(), last.begin(), last.end());
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
Intor Dly_to_Mly_DateList(Intor Dly_DateList)
{
    Intor Mly_DateList(0);
    for(int i=0; i<Dly_DateList.size(); i++)
        if (Dly_DateList[i] - Dly_DateList[i-1] > 40)
            Mly_DateList.push_back(Dly_DateList[i-1]);
    Mly_DateList.push_back(Dly_DateList[Dly_DateList.size()-1]);
    return Mly_DateList;
}
Matrix combine_First_with_Mly_MC(Matrix MC, Matrix fMC)
{
    Matrix combined_MC;
    int i=0, j=0, n_MC=MC.size(), n_fMC=fMC.size();
    int last_ID=9999;

    while(i < n_MC-1 || j < n_fMC-1) //todo: <= or <
    {
        while  (MC[i][0] <= last_ID && i<n_MC -1)   i++;
        while (fMC[j][0] <= last_ID && j<n_fMC-1)  j++;

        if    (MC[i][0] <= fMC[j][0] && MC[i][0] > last_ID)
        {
            combined_MC.push_back(MC[i]);
            last_ID = MC[i][0];
        }
        else if (fMC[j][0] > last_ID)
        {
            Vector stock(3);
            for (int k = 0; k<3; ++k)   stock[k] = fMC[j][k];
            combined_MC.push_back(stock);
            last_ID = stock[0];
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
        if(i == Rs_Dates.size()-2)
            cout << "i-4";
        if(pMC2[j][0] == Rs_Dates[i][0])
        {
            pMC_comb.push_back(pMC2[j]);
            j++;
        }
        else {
            cout << "ayy";
            pMC_comb.push_back({Rs_Dates[i][0]*1.0, Rs_Dates[i][1]/100.0, 0});
        }
    }
    return pMC_comb;
}
void defineFilePaths(string& incr, string& Exo_FilePath, string& Proccessed_FilePath, string& Proccessed_FilePath_inc, string& Proccessed_Dly, string& Proccessed_Mly, string& Proccessed_Yly)
{
    Exo_FilePath = "Data/Input/Exo_Files/";
    Proccessed_FilePath = "Data/Input/Processed_Files/";
    Proccessed_Dly = Proccessed_FilePath+"Dly/";
    Proccessed_Mly = Proccessed_FilePath+"Mly/";
    Proccessed_Yly = Proccessed_FilePath+"Yly/";
    Proccessed_FilePath_inc = Proccessed_FilePath + incr + "/";
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
void remove_substring(string& str, string sub_str)
{
    std::string::size_type i = str.find(sub_str);

    if (i != std::string::npos)
        str.erase(i, sub_str.length());
}
void remove_substring(vector<string>& v_str, string sub_str)
{
    for(auto& str:v_str)
        remove_substring(str, sub_str);
}
Vector int_to_double(Intor intVec)
{
    vector<double> doubleVec(intVec.begin(), intVec.end());
    return doubleVec;
}