#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>        //min function

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
Matrix Beta_Alpha_Calculate(Matrix DR, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, const Intor& iPeriod, int minTradingDaysBeforePeriod);
void TestCalculations(int max, double max_ratio, int minTradingDays);
void RunCalculations(int max, double max_ratio, int minTradingDays);

int main()
{
    //Create_Files_from_DR(-1);
    //TestCalculations(-1, 0.4, 75);
    RunCalculations(100, 0.4, 75);
}

void TestCalculations(int max, double max_ratio, int minTradingDays)
{
    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,max_ratio, minTradingDays);

    //Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    Intrix iPeriod = Load_Intrix("Yearly_iPeriods.txt", -1);

    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    cout << "Files was Loaded.\n\n";

    Vector beta(0), alpha(0), akk_return(0), PERMNO(0);
    //vector<Matrix> Data;

    for (int i = 0; i < iPeriod.size(); ++i)
    {
        //Calculate {beta, alpha, return, ID}
        Matrix beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod[i], minTradingDays);

        //Store data
        beta.insert(beta.end(), beta_alpha_return[0].begin(), beta_alpha_return[0].end());
        alpha.insert(alpha.end(), beta_alpha_return[1].begin(), beta_alpha_return[1].end());
        akk_return.insert(akk_return.end(), beta_alpha_return[2].begin(), beta_alpha_return[2].end());
        PERMNO.insert(PERMNO.end(), beta_alpha_return[3].begin(), beta_alpha_return[3].end());
        //Data.push_back(beta_alpha_return);
    }
    //Print means
    cout << "mean(beta) = " << Vector_Average(beta) << endl;
    cout << "mean(alpha) = " << Vector_Average(alpha) << endl;
    cout << "mean(akk_return) = " << Vector_Average(akk_return) << endl;

    //Save Data
    Save_Vector("Data_beta.txt",beta);
    Save_Vector("Data_alpha.txt",alpha);
    Save_Vector("Data_akk_return.txt",akk_return);
    Save_Vector("Data_PERMNO.txt",PERMNO);
}
void RunCalculations(int max, double max_ratio, int minTradingDays)
{
    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Intrix iPeriods = Load_Intrix("iPeriods.txt", -1);
    Intrix iPeriods = Load_Intrix("Monthly_iPeriods.txt", -1);

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector("sp500.txt");
    Vector riskFree = Load_Vector("DailyRiskFreeReturn.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");

    cout << "Files was Loaded.\n\n";

    Vector beta(0), alpha(0), akk_return(0), PERMNO(0);
    Matrix beta_alpha_return;
    for (auto & iPeriod : iPeriods)
    {
        //Calculate Data
        beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, 130);
        //Store Data
        beta.insert(beta.end(), beta_alpha_return[0].begin(), beta_alpha_return[0].end());
        alpha.insert(alpha.end(), beta_alpha_return[1].begin(), beta_alpha_return[1].end());
        akk_return.insert(akk_return.end(), beta_alpha_return[2].begin(), beta_alpha_return[2].end());
        PERMNO.insert(PERMNO.end(), beta_alpha_return[3].begin(), beta_alpha_return[3].end());
    }
    //Save Data
    Save_Vector("DataRun_beta.txt",beta);
    Save_Vector("DataRun_alpha.txt",alpha);
    Save_Vector("DataRun_akk_return.txt",akk_return);
    Save_Vector("DataRun_PERMNO.txt",PERMNO);
}
Matrix Beta_Alpha_Calculate(Matrix DR, Intrix iDates, const Vector& sp500, const Vector& riskFree, const Intor& Dates, const Intor& iPeriod, int minTradingDaysBeforePeriod)
{
    Vector cov(0), sp500_Var(0), beta(0), stock_return_akk(0), sp500_return_akk(0), alpha(0), PERMNO(0);
    int iStart_sp500, iEnd_sp500, iStart_DR, iLength, Active_days;
    double stockBeta, stockAlpha, sp500_Return_akk, stock_Return_akk, RiskFree_Return_akk, Active_years;
    Vector akk_stock_sp500_RiskFree_Return;
    int emptyCount;

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
        emptyCount=0;   if(iLength > 0) for(int j = iStart_DR; j < iStart_DR+iLength; ++j)  if(DR[i][j] == -2)  emptyCount++;
        Active_days = iLength - emptyCount;
        Active_years = (DaysBetween(Dates[iStart_sp500], Dates[iStart_sp500+Active_days]) + 0.0) / 365.24;

        //Yderligere krav om at aktien ikke dør indenfor meget få dage //TODO: skaber bias, skal helst fjernes/mindskes på sigt
        if((Active_days) < 10) continue;

        //Calculate beta
        stockBeta = Calculate_stockBeta(DR[i], sp500, iStart_sp500, iEnd_sp500, iStart_DR);

        //Calculate akk. returns
        akk_stock_sp500_RiskFree_Return = Calculate_akk_Return(DR[i], sp500, riskFree, iLength, iStart_DR, iStart_sp500);
        stock_Return_akk = akk_stock_sp500_RiskFree_Return[0];
        sp500_Return_akk = akk_stock_sp500_RiskFree_Return[1];
        RiskFree_Return_akk = akk_stock_sp500_RiskFree_Return[2];

        //Calculate alpha
        stockAlpha = Calculate_StockAlpha(stockBeta, stock_Return_akk, sp500_Return_akk, RiskFree_Return_akk, Active_years);

        //Save data to vectors
        beta.push_back(stockBeta);
        alpha.push_back(stockAlpha);
        stock_return_akk.push_back(stock_Return_akk);
        sp500_return_akk.push_back(sp500_Return_akk);
        PERMNO.push_back(DR[i][0]);
    }
    return {beta, alpha, stock_return_akk, PERMNO};
}