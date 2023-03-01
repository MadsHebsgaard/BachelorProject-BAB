#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once
#define BACHELOR_CALCULATIONS_H

double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR);      //TODO
double Calculate_stockReturn(Vector stock, int iLength, int iStart_DR);             //TODO
double CovSP500(Vector Stock, Vector sp500, int iStart, int iEnd, int iStart_DR);   //TODO: Fix start date and end date
double Var(Vector sp500, int Start_number, int End_number);                         //TODO: Fix start date and end date
int day_after1926(int date);
Intor Column_of_Intrix(Intrix A, int n);
int Find_date_integer(int date, Intor& Date_list);
double Vector_Average(const Vector& v);
Vector Matrix_Column(Matrix A, int j);


double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR)
{
    double stockCov = CovSP500(stock, sp500, iStart, iEnd, iStart_DR);
    double stockSP500var = Var(sp500, iStart, iEnd);
    return (stockCov/stockSP500var);
}
double Calculate_stockReturn(Vector stock, int iLength, int iStart_DR)
{
    double Return=1;
    int trading_Days = 0;
    for (int i = iStart_DR; i < iStart_DR + iLength; ++i)
    {
        if(stock[i] == -2) continue;
        trading_Days++;
        Return*=(1 + stock[i]);
    }
    return Return - 1;                                    //TODO: Return type
    //return pow(Return,263.5104/trading_Days)-1;        //
}
double CovSP500(Vector Stock, Vector sp500, int iStart, int iEnd, int iStart_DR)
{
    int d_iDates = iEnd - iStart + 1;
    double Mean_Stock = 0, Mean_sp500 = 0, Covv = 0;
    int trading_Days = 0;

    //Mean
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i] != -2)
        {
            trading_Days++;
            Mean_Stock += Stock[i + iStart_DR];
            Mean_sp500 += sp500[i + iStart];
        }
    }
    Mean_Stock /= trading_Days;
    Mean_sp500 /= trading_Days;

    //CovSP500
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i] != -2)
            Covv += (Stock[i+iStart_DR] - Mean_Stock) * (sp500[i + iStart] - Mean_sp500);
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
    int træk, zero_total;
    int max_træk = 6;
    int replace = -2;
    int continue_amount = 0;

    for (auto & stock : A)
    {
        zero_total = 0;
        træk = 0;
        if (stock.size() < minTradingDays)   {continue_amount++; continue;}
        for (int j = 0; j < stock.size(); ++j)
        {
            if(stock[j] == 0 or stock[j] == replace)    {træk++; zero_total++;}
            else                 træk = 0;

            if      (træk == max_træk)  for(int k = 1; k <= max_træk; ++k) stock[j - max_træk + k] = replace;
            else if (træk > max_træk) stock[j] = replace;
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
