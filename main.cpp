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

void Create_all_CSV_files(string incr)
{
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    double max_ratio = 0.4;
    int minTradingDays = 10; //10 / 140
    int Rs_size = -1;   // load all
    Matrix Rs = Load_Rs_Compressed(Proccessed_FilePath_incr + "Rs.txt", Rs_size);
    Matrix Rs_i = Edit_Rs(Rs, max_ratio);

    string filename = "run";
    bool logarithm = false;

    for (int i = 0; i<2; ++i) {
        Simple_run(incr, filename, minTradingDays, Rs, logarithm);
        PrePost_run(incr, filename, Rs, logarithm);

        vector<string> calc_methods = {"RF", "DT", "RB"};
        for (auto& str:calc_methods) {
            BackTest_run(incr, "x_"+str+"_" + filename, filename, Rs, {-2, 1, 4}, str, false);
            BackTest_run(incr, "x_"+str+"_" + filename + "_w_m", filename, Rs, {-2, 4}, str, true);
            BackTest_run(incr, "x_"+str+"_" + filename + "_w", filename, Rs, {-2, 1, 4}, str, true);
            BackTest_run(incr, "x_"+str+"_" + filename +"_m", filename, Rs, {-2, 4}, str, false);
            BackTest_exotic(incr, "x_"+str+"_" + filename + "_momentum", filename, Rs, {-0.4, 0.5, 1}, str, "momentum", false);
            BackTest_exotic(incr, "x_"+str+"_" + filename + "_MC", filename, Rs, {0, 0.5, 20}, str, "MC", false);
        }
        filename = filename+"_log";
        logarithm=true;
    }
}

void create_core_CSV_files(string incr)
{
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    double max_ratio = 0.4;
    int minTradingDays = 10; //10 / 140
    int Rs_size = -1;   // load all
    Matrix Rs = Load_Rs_Compressed(Proccessed_FilePath_incr + "Rs.txt", Rs_size);
    Matrix Rs_i = Edit_Rs(Rs, max_ratio);

    string filename = "run", method="RF";
    bool logarithm = true;

    Simple_run(incr, filename, minTradingDays, Rs, logarithm);
    PrePost_run(incr, filename, Rs, logarithm);
    BackTest_run(incr, filename, filename, Rs, {-2, 1, 4}, method, logarithm);
    BackTest_run(incr, filename + "_w_m", filename, Rs, {-2, 4}, method, logarithm);
}
void play(string incr)
{
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    double max_ratio = 0.4;
    int minTradingDays;
    if(incr=="Mly") minTradingDays=10;
    else            minTradingDays=140;

    int Rs_size = -1;   // load all
    Matrix Rs = Load_Rs_Compressed(Proccessed_FilePath_incr + "Rs.txt", Rs_size);
    Matrix Rs_i = Edit_Rs(Rs, max_ratio);

    string filename = "run", method="RF";
    //bool logarithm = true;

    BackTest_exotic(incr, "qq_" + filename + "_score_w_m", filename, Rs, {-2, 4}, method, "score", true);

    //Simple_run(incr, filename, minTradingDays, Rs, logarithm);
    //PrePost_run(incr, filename, Rs, logarithm);
    //BackTest_run(incr, filename, filename, Rs, {-2, 1, 4}, method, logarithm);
    //BackTest_run(incr, filename + "_w_m", filename, Rs, {-2, 4}, method, logarithm);
}
int main()
{
    //Setup_all_files();  //Sets up all files. Turn off when you have done this sucessfully once.
    //Create_all_CSV_files("Mly");   //Mly  Dly (Creates all files for monthly / daily)
    //create_core_CSV_files("Mly");   //Creates all core monthly BAB-files
    //create_core_CSV_files("Dly");   //Creates all core daily BAB-files
    play("Dly");
    return 0;
}