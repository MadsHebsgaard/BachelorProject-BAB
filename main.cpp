#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "load.h"           //Load functions
#include "save.h"           //Save functions
#include "ui.h"             //UI stuff
#include "calculations.h"   //any calculation
#include "testing.h"        //to be used for testing only

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

//to be changed
Matrix Beta_Alpha_Calculate(Matrix DR, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& iPeriod, int minTradingDaysInPeriod);
void TestCalculations(int max);
void RunCalculations(int max);

//hva så røvhuller

int main()
{
    //sp500_Dates_to_monthly_dDate_periods();
    //Create_Files_from_DR(-1);
    TestCalculations(-1);
    //RunCalculations(500);
    return 0;
}

void TestCalculations(int max)
{
    //max = 500;
    int minTradingDays = 75;

    //Matrix DR_ny = Load_DR_Compressed("DR_Compressed.txt", max);
    //Intrix StockDays = Load_StockDays_Compressed("DR_StockDays.txt", true, max);
    //Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);
    //Intrix iPeriod = Load_Intrix("Yearly_iPeriods.txt", -1);

    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,0.3,minTradingDays);

    Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);
    Vector beta(0); //TODO: One for mu and alpha too

    //Beta
    int maxPeriods = iPeriod.size();
    maxPeriods = 1;
    for (int i = 0; i < maxPeriods; ++i)
    {
        Matrix beta_mean = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, iPeriod[i], 130);
        Data.push_back(beta_mean);

        if(i < 15)
        {
            cout << "n = "     << setw(5) << beta_mean[0].size();
            cout << ",  Avg_B = " << setw(5) << Vector_Average(beta_mean[0]);
            cout << ",  Avg_a = " << setw(5) << Vector_Average(beta_mean[1]) << endl;
        }
        //beta.insert(beta.end(), beta_mean[0].begin(), beta_mean[0].end());
    }
    Save_Vector("Data_beta.txt",Data[0][0]);
    Save_Vector("Data_alpha.txt",Data[0][1]);
}
void RunCalculations(int max)
{
    //max = 5000;
    int minTradingDays = 75;

    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,0.3,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);

    //Load sp500 stock returns
    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);

    //Beta
    for (int i = 0; i < iPeriod.size(); ++i)
    {
        Matrix beta_mean = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, iPeriod[i], 8);
        Data.push_back(beta_mean);

        cout << setw(9) << "n = " << beta_mean[0].size();
        cout << setw(18) << "Avg_B = " << Vector_Average(beta_mean[0]) << endl;
    }
}
Matrix Beta_Alpha_Calculate(Matrix DR, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& iPeriod, int minTradingDaysInPeriod)
{
    //TODO: (done) skip the first minTradingDaysInPeriod of the stock
    Vector cov(0), sp500_Var(0), beta(0), stock_return_akk(0), sp500_return_akk(0), alpha(0);
    int iStart_sp500, iEnd, iStart_DR, iLength, Active_days;
    double stockBeta, stockAlpha, sp500_Return_akk, sp500_avg_Return, riskFree_avg_Return, stock_Return_akk, RiskFree_Return_akk;
    Vector akk_stock_sp500_RiskFree_Return;
    int emptyCount;

    for (int i = 0; i < DR.size(); ++i)
    {
        //Determen index start and end in sp500 and index to start in DR[i]
        iStart_DR = 1;

        if(iDates[i][1] < iPeriod[0])               iStart_DR = 1 + (iPeriod[0] - iDates[i][1]);    //index where period start in DR[i]
        if(iStart_DR < 1+minTradingDaysInPeriod)    iStart_DR = 1 + minTradingDaysInPeriod;       //index where period min start in DR[i]

        iStart_sp500 = iDates[i][1] + iStart_DR-1;      //index where sp500 start
        iEnd = min(iDates[i][2], iPeriod[1]);               //sp500 date to end
        iLength = iEnd - iStart_sp500 + 1;                        //iDates total

        //Check if stock is relevant
        emptyCount=0;
        if(iLength > 0) for (int j = iStart_DR; j < iStart_DR+iLength; ++j)     if(DR[i][j] == -2)    emptyCount++;
        Active_days = iLength - emptyCount;
        if((Active_days) < 10) continue;    //TODO: Burde mÃ¥ske nÃ¦rmere vÃ¦re <2 / <3 / <4

        //Calculate beta
        stockBeta = Calculate_stockBeta(DR[i], sp500, iStart_sp500, iEnd, iStart_DR);

        //Calculate return
        akk_stock_sp500_RiskFree_Return = Calculate_akk_Return(DR[i], sp500, riskFree, iLength, iStart_DR, iStart_sp500);
        stock_Return_akk = akk_stock_sp500_RiskFree_Return[0];
        sp500_Return_akk = akk_stock_sp500_RiskFree_Return[1];
        RiskFree_Return_akk = akk_stock_sp500_RiskFree_Return[2];

        //Non akk. returns
        sp500_avg_Return = pow(1 + sp500_Return_akk, 263.5104 / Active_days) - 1;
        riskFree_avg_Return = pow(1 + RiskFree_Return_akk, 263.5104 / Active_days) - 1;

        stockAlpha = stock_Return_akk - (pow(1 + stockBeta * (sp500_avg_Return - riskFree_avg_Return) + riskFree_avg_Return, Active_days / 263.5104) - 1); //TODO: /T

        //cout << " " << riskFree_avg_Return << " ";
        //cout << iDates[i][0] << setw(5) << Active_days << ", sp500_Return_akk =" << setw(12) << sp500_Return_akk;
        //cout << ", stock_Return_akk =" << setw(12) << stock_Return_akk;
        //cout << ", stockAlpha =" << setw(12) << stockAlpha << endl;

        beta.push_back(stockBeta);
        alpha.push_back(stockAlpha);
        stock_return_akk.push_back(stock_Return_akk);
        sp500_return_akk.push_back(sp500_Return_akk);
    }

    /*
    cout << endl;
    int min_stock_return_akk_index = min_Val_Index(stock_return_akk);
    int min_sp500_return_akk_index = min_Val_Index(sp500_return_akk);
    cout << min_stock_return_akk_index << ",  min stock_return_akk = " << stock_return_akk[min_stock_return_akk_index] << endl;
    cout << min_sp500_return_akk_index << ",  min sp500_return_akk = " << sp500_return_akk[min_sp500_return_akk_index] << endl;
    */

    return {beta, alpha};
}
