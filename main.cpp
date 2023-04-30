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
#include "backtest.h"       //Backtesting BAB with 'FindPortfolios'
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
using Tensor5 = vector<Tensor4>;

int main()
{
    //Setup_all_files();

    string incr = "Mly";  //Mly  Dly
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);
    //Process_Files(incr, false, false);

    Matrix Rs = Load_Rs_Compressed(Proccessed_FilePath_incr + "Rs.txt", -1);

    Simple_Calculations(incr, "Run", 0.3, 12, Rs);
    PrePost_Calculations(incr, "Run", 0.3, Rs);
    FindPortfolios(incr, "Run", "Run", Rs);

    return 0;
}
