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
#include "createData.h"     //Create era/period data
#include "load_output.h"    //Load output files/data

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;

//Process_Files();
//Era_Calculations("test", 0.3, 100, 1, DR);
//Period_Calculations("test", 0.3, 100, DR, 1);
//Era_PrePost_Calculations("test", 0.4, 100, 1, DR);
//Era_Period_PrePost_Calculations("test", 0.4, 100, 1, DR);


//Tensor5 Data = Load_Era_PrePost_Period("run");
//Tensor4 Data = Load_Era_n_PrePost_Period("run",1);
//vector<Matrix> Data = Load_Era_PrePost("Test", 1);



int main()
{
    //Process_Files();
    //Matrix DR = Load_DR_Compressed("Data/Input/Processed_Files/DR_Compressed.txt", -1);
    //Era_Period_PrePost_Calculations("run", 0.3, 100, 1, DR);

    //Era_Calculations("run_test60", 0.3, 100, 1, DR);
    //Period_Calculations("run_test60", 0.3, 100, DR, 1);
    //Era_PrePost_Calculations("runtest", 0.3, 100, 5, DR);
    //Era_PrePost_Period_Calculations("run", 0.3, 100, 5, DR);
    //Tensor5 Data = Load_Era_PrePost_Period("run");


    Tensor5 Data = Load_Era_Period_PrePost("run");
    BackTesting(Data, 1);
}
