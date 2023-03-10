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

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

//to be changed
void synchronousCalculations(string folderName, double max_ratio, int minTradingDays, int n_periods, Matrix DR);
void DuoPeriod_Calculations(string folderName, double max_ratio, int minTradingDays, int n_periods, Matrix DR);

//TODO: Needed files list
//TODO: Easily create rest of files

int main()
{
    //Process_Files();
    Matrix DR = Load_DR_Compressed("Data/Input/Processed_Files/DR_Compressed.txt", 50);
    synchronousCalculations("Mads", 0.4, 75, 6, DR);
    DuoPeriod_Calculations("Maria", 0.4, 75, 6, DR);
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

void DuoPeriod_Calculations(string folderName, double max_ratio, int minTradingDays, int n_periods, Matrix DR)
{
    string Exo_FilePath = "Data/Input/Exo_Files/";
    string Proccessed_FilePath = "Data/Input/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    //Load iPeriods and create Era_List
    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    vector<Intrix> Era_List = SplitPeriods(iPeriods, n_periods, true);

    cout << "DuoPeriod_Calculations: Files was Loaded for \"" << folderName << "\".\n\n";
    mkdir("Data/Output");
    mkdir("Data/Output/Double");
    folderName = "Data/Output/Double/" + folderName;
    mkdir(folderName.c_str());
    vector<Matrix> twoPeriod_Data;
    vector<vector<Matrix>> tP_DataSet(Era_List.size(),vector<Matrix>(2,Matrix(0,Vector(0))));

    for (int Era = 0; Era < Era_List.size(); ++Era)
    {
        Matrix beta(2,Vector(0)), alpha(2,Vector(0)), PERMNO(2,Vector(0)), akk_return(2,Vector(0)), akk_sp500(2,Vector(0)), akk_riskFree(2,Vector(0));
        for (int Period = 0; Period < Era_List[Era].size()-1; Period++)
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            twoPeriod_Data = TwoPeriod_Calc(DR_ny, iDates, sp500, riskFree, Dates, Era_List[Era][Period], Era_List[Era][Period + 1]);
            push_back(tP_DataSet[Era], twoPeriod_Data);
        }
        //Create directory and path for era
        string dirName = folderName + "/Era" + "_" + to_string(Era + 1);
        mkdir(dirName.c_str());
        string preDirName = dirName + "/Pre_Period";
        string periodDirName = dirName + "/Period";

        //Create files in directory for the era
        vector<string> paths {preDirName, periodDirName};
        vector<string> fileNames = {"/beta.txt", "/alpha.txt", "/akk_return.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
        for (int i = 0; i < paths.size(); ++i) {
            mkdir(paths[i].c_str());
            for (int file = 0; file < 6; ++file)    Save_Vector(paths[i] + fileNames[file],tP_DataSet[Era][i][file]);
        }
    }
}
