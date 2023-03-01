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
Matrix Alpha_Beta_Mu(Matrix DR, Intrix iDates, Vector sp500, const Intor& iPeriod, int minTradingDaysInPeriod);
void TestCalculations(int max);
void RunCalculations(int max);

int main()
{
    //sp500_Dates_to_monthly_dDate_periods();
    //Create_Files_from_DR(-1);
    TestCalculations(5000);
    //RunCalculations(500);
    return 0;
}

void TestCalculations(int max)
{
    //max = 500;
    int minTradingDays = 75;

        //DR with condition for inclusion
    //Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    //Matrix DR_ny = Edit_DR(DR,0.3,minTradingDays);

    Matrix DR_ny = Load_DR_Compressed("DR_Compressed.txt", max);
    //Intrix StockDays = Load_StockDays_Compressed("DR_StockDays.txt", true, max);
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    //Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);
    Intrix iPeriod = Load_Intrix("Yearly_iPeriods.txt", -1);
    Vector sp500 = Load_Vector("sp500.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");


    //iDates with same stocks as DR_ny
    //Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    //iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

        //Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    //Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);

    Vector beta(0); //TODO: One for mu and alpha too

    //Beta
    for (int i = 0; i < iPeriod.size(); ++i)
    {
        Matrix beta_mean = Alpha_Beta_Mu(DR_ny, iDates, sp500, iPeriod[i], 130);
        Data.push_back(beta_mean);

        if(i < 15)
        {
            cout << setw(9) << "n = " << beta_mean[0].size();
            cout << setw(18) << "Avg_B = " << Vector_Average(beta_mean[0]) << endl;
        }
        beta.insert(beta.end(), beta_mean[0].begin(), beta_mean[0].end());
    }
    //cout << "mean beta = " << Vector_Average(Data[0]) << endl;
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

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);

    //Beta
    for (int i = 0; i < iPeriod.size(); ++i)
    {
        Matrix beta_mean = Alpha_Beta_Mu(DR_ny, iDates, sp500, iPeriod[i], 8);
        Data.push_back(beta_mean);

        cout << setw(9) << "n = " << beta_mean[0].size();
        cout << setw(18) << "Avg_B = " << Vector_Average(beta_mean[0]) << endl;
    }
}
Matrix Alpha_Beta_Mu(Matrix DR, Intrix iDates, Vector sp500, const Intor& iPeriod, int minTradingDaysInPeriod)
{
    Vector cov(0);
    Vector sp500_Var(0);
    Vector beta(0);
    Vector returns(0);
    Vector alpha(0);
    int iStart, iEnd, iStart_DR, iLength;
    double stockBeta, stockReturn, stockAlpha;
    int emptyCount;
    //double stockCov, stockSP500var;



    for (int i = 0; i < DR.size(); ++i)
    {
        //Determen index start and end in sp500 and index to start in DR[i]
        if(iDates[i][1] < iPeriod[0])
        {
            iStart = iPeriod[0];                            //sp500 date to start
            iStart_DR = 1 + (iPeriod[0] - iDates[i][1]);    //index where period start in DR[i]
        }
        else
        {
            iStart = iDates[i][1];
            iStart_DR = 1;
        }
        iEnd = min(iDates[i][2], iPeriod[1]);               //sp500 date to end
        iLength = iEnd - iStart + 1;                        //iDates total

        //Check if stock is relevant
        emptyCount=0;
        if(iLength > 0) for (int j = iStart_DR; j < iStart_DR+iLength; ++j)     if(DR[i][j] == -2)    emptyCount++;
        if((iLength - emptyCount) < minTradingDaysInPeriod) continue;

        //Calculate beta
        stockBeta = Calculate_stockBeta(DR[i], sp500, iStart, iEnd, iStart_DR);
        beta.push_back(stockBeta);

        //Calculate return
        stockReturn = Calculate_stockReturn(DR[i], iLength, iStart_DR);
        returns.push_back(stockReturn);

        double rf_return = 0.02;    //TODO: insert periods risk free return as an argument in the function Alpha_Beta_Mu
        stockAlpha = (stockReturn-rf_return)/stockBeta;
        alpha.push_back(stockAlpha);

        //cout << "stockReturn = " << setw(10) << stockReturn << " i = " << i << endl;
    }


    double sp500_mean = Calculate_stockReturn(sp500, (iPeriod[1] - iPeriod[0]) + 1, iPeriod[0]);
    //cout << setw(18) << "sp_r = " << sp500_mean;
    //cout << "\nreturns = \n";    Print_Vector(returns);
    //cout << setw(18) << "avg_r = " << Vector_Average(returns);
    //if(DR.size() < 50) {cout << "cov = \n"; Print_Vector(cov); cout << "sp500_Var = \n"; Print_Vector(sp500_Var);}
    return {beta, returns, alpha};
}
