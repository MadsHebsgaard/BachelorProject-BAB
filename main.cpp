//C++ standard functions
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>        //min function
#include <sys/stat.h>       // For mkdir()

//Our functions
#include "process.h"        //Data creation functions
#include "createData.h"     //Create 'Simple' and 'PrePost' data
#include "backtest.h"       //Backtesting BAB with 'BackTest_run'
#include "load.h"           //Load functions
#include "save.h"           //Save functions
#include "calculations.h"   //any calculation
#include "load_output.h"    //Load output files/data
#include "ui.h"             //UI/console functions
#include "testing.h"        //to be used for testing only
#include "LegacyCode.h"     //Old code no longer in use

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;

//todo:
//(proccesed files) Add yly RFR
//(proccesed files) Add yly sp500
//(backtesting)     Antal aktier i hver portefølge start og slut
//(Setup_all_files) Fix mly iPeriods
//Check for DR

int main()
{

    string incr = "Mly";  //Mly  Dly
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    //Setup_all_files();

    double max_ratio = 0.4;
    int minTradingDays = 10;

    Matrix Rs = Load_Rs_Compressed(Proccessed_FilePath_incr + "Rs.txt", -1);
    Rs = Edit_DR(Rs, max_ratio, minTradingDays);    //Dly problem when TD is low and R = 1. May cause some bias, so turned off

    string filename = "run2";
    Simple_run(incr, filename, minTradingDays, Rs);
    PrePost_run(incr, filename, Rs);

    BackTest_run(incr, "RF_"+ filename, filename, Rs, {-0.5, 1, 1, 3}, "RF");
    BackTest_run(incr, "DT_"+ filename, filename, Rs, {-0.5, 1, 1, 3}, "DT");

    return 0;
}