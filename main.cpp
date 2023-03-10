#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>        //min function
#include <sys/stat.h>       // For mkdir()

#include "load.h"           //Load functions
#include "save.h"           //Save functions
#include "ui.h"             //UI stuff
#include "calculations.h"   //any calculation
#include "testing.h"        //to be used for testing only
#include "createData.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;


int main()
{
    //Process_Files();
    Matrix DR = Load_DR_Compressed("Data/Input/Processed_Files/DR_Compressed.txt", 100);
    SingleEra_Calculations("Run", 0.3, 100, 5, DR);
    //DuoPeriod_Calculations("Run", 0.3, 100, 5, DR);
    //SinglePeriod_Calculations("Run", 0.3, 100, DR, 1);



}