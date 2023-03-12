#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "load.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once


void push_back(vector<Matrix>& First, const vector<Matrix>& Second);
void push_back(Matrix& First, const Matrix& Second);

double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR);      //TODO
pair<double,double> Calculate_Return(Vector stock, int iLength, int iStart_DR);             //TODO
double CovSP500(Vector Stock, Vector sp500, int iStart_sp500, int iEnd_sp500, int iStart_DR);   //TODO: Fix start date and end date
double Var(Vector sp500, int Start_number, int End_number);                         //TODO: Fix start date and end date
int day_after1926(int date);
Intor Column_of_Intrix(Intrix A, int n);
int Find_date_integer(int date, Intor& Date_list);
double Vector_Average(const Vector& v);
Vector Matrix_Column(Matrix A, int j);
void Save_Vector(const string& fn, const Vector& v);
Intor Load_Intor(const string& fn);

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
int Find_date_integer(int date, Intor& Date_list)
{
    int index_apprx = day_after1926(date);
    while(Date_list[index_apprx] < date)    index_apprx += 15;
    while(Date_list[index_apprx] > date)    index_apprx -= 1;
    return index_apprx;
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
Intrix Dates_to_iDates(Intrix Dates, Intor True_Dates, int start_j)
{
    for (auto & Date : Dates)
    {
        for (int j = start_j; j < Date.size(); ++j)
        {
            Date[j] = Find_date_integer(Date[j], True_Dates);
        }
    }
    return Dates;
}
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
        if((zero_total+0.0) / stock.size() < max_ratio)  {DR_ny.push_back(stock);           /*cout << endl << "zero_total:" << zero_total << "  size: " << A[stock].size();*/}
        //else cout << endl << stock << "   " << A[stock][0] << "  " << (zero_total+0.0)/A[stock].size();
    }

    cout << "\ncontinue_amount = " << continue_amount;
    cout << "\nDR_ny size = " << DR_ny.size();
    cout << "\nZero ratio too high = " << A.size() - DR_ny.size() - continue_amount;
    cout << "\n\n";

    return DR_ny;
}
Intrix Remove_Missing_ID_Intrix(Intrix A, Vector v)
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
Intrix dPeriods_From_iPeriod(Intrix iPeriod)
{
    Intor Dates = Load_Intor("DateList.txt");
    for(int i =0; i<iPeriod.size(); i++)
    {
        iPeriod[i][0] = Dates[iPeriod[i][0]];
        iPeriod[i][1] = Dates[iPeriod[i][1]];
    }
    return iPeriod;
}
Intrix Yearly_iPeriods(Intor DateList)
{
    int i = 0, periodStart = DateList[0];

    i = 0;
    periodStart = DateList[0];
    Intrix YearlyPeriods(0);
    while (DateList[i] < 99999999)
    {
        if((DateList[i] - periodStart) > 4000) //40 for monthly
        {
            YearlyPeriods.push_back({periodStart, DateList[i - 1]});
            periodStart = DateList[i];
        }
        i++;
    }
    Intrix Yearly_iPeriods = Dates_to_iDates(YearlyPeriods, DateList, 0);
    return Yearly_iPeriods;
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
Matrix Beta_Alpha_Calculate(Matrix DR, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, const Intor& iPeriod, int minTradingDaysBeforePeriod)
{
    Vector beta(0), stock_return_akk(0), sp500_return_akk(0), alpha(0), PERMNO(0), riskFree_Return_akk(0);
    int iStart_sp500, iEnd_sp500, iStart_DR, iLength, Active_days;
    double stockBeta, stockAlpha, sp500_Return_akk, stock_Return_akk, RiskFree_Return_akk, Active_years;
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
            if(DR[i][j] == -2)
                emptyCount++;
        Active_days = iLength - emptyCount;
        if(Active_days < 1) continue;

        //Active_TradingDaysBeforePeriod
        Active_TradingDaysBeforePeriod=0;
        for(int j = iStart_DR-minTradingDaysBeforePeriod; j < iStart_DR; ++j)
            if(DR[i][j] > -1.5 && DR[i][j] != 0)
                Active_TradingDaysBeforePeriod++;

        //Yderligere krav om at aktien er blevet handlet før tid
        if(1+Active_TradingDaysBeforePeriod < minTradingDaysBeforePeriod/3.0) continue;

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

        //Save data to vectors
        beta.push_back(stockBeta);
        alpha.push_back(stockAlpha);
        stock_return_akk.push_back(stock_Return_akk);
        sp500_return_akk.push_back(sp500_Return_akk);
        riskFree_Return_akk.push_back(RiskFree_Return_akk);
        PERMNO.push_back(DR[i][0]);
    }
    //double avg_sp500 = Calculate_akk_Return(sp500, sp500, riskFree, riskFree.size(), 0, 0)[1];
    //double avg_riskFree = Calculate_akk_Return(sp500, sp500, riskFree, riskFree.size(), 0, 0)[2];
    return {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk};
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
    int j = 0;
    for (int i = 0; i < A[0].size(); ++i)
    {
        while(A[ID_row][i] > B[ID_row][j])
        {
            j++;
            if(j>=B.size()) break;
        }
        if(A[ID_row][i] == B[ID_row][j])
        {
            A_Keep.push_back(i);
            B_Keep.push_back(j);
            j++;
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
vector<Matrix> TwoPeriod_Calc(Matrix DR, Intrix iDates, Vector sp500, Vector riskFree, Intor Dates, Intor Pre_iPeriod, Intor iPeriod)
{
    vector<Matrix> A(0);
    Matrix Pre_Period_values = Beta_Alpha_Calculate(DR, iDates, sp500, riskFree, Dates, Pre_iPeriod, 0);
    Matrix Period_values = Beta_Alpha_Calculate(DR, iDates, sp500, riskFree, Dates, iPeriod, 20); //20 is somewhat arbitrary, should make beta less likely to be numerically huge
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