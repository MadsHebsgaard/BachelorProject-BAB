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
void RunCalculations(string folderName, double max_ratio, int minTradingDays, int n_periods, Matrix DR);

//TODO: Needed files list
//TODO: Easily create rest of files

int main()
{
    //HowToGetStarted();
    //Create_Files();
    Matrix DR = Load_DR_Compressed("Input_Data/Processed_Files/DR_Compressed.txt", 500);
    RunCalculations("Run", 0.4, 75, 95, DR);
}

void RunCalculations(string folderName, double max_ratio, int minTradingDays, int n_periods, Matrix DR)
{
    string Exo_FilePath = "Input_Data/Exo_Files/";
    string Proccessed_FilePath = "Input_Data/Processed_Files/";

    //DR with condition for inclusion
    Matrix DR_ny = Edit_DR(DR,max_ratio,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix(Proccessed_FilePath+"DR_iDates.txt", -1);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Load sp500 and riskFree returns & SP500 dates
    Vector sp500 = Load_Vector(Exo_FilePath+"sp500.txt");
    Vector riskFree = Load_Vector(Proccessed_FilePath+"riskFreeReturn.txt");
    Intor Dates = Load_Intor(Exo_FilePath+"DateList.txt");

    cout << "Files was Loaded.\n\n";

    Matrix beta(n_periods,Vector(0)), alpha(n_periods,Vector(0)), PERMNO(n_periods,Vector(0));
    Matrix akk_return(n_periods,Vector(0)), akk_sp500(n_periods,Vector(0)), akk_riskFree(n_periods,Vector(0));
    Matrix beta_alpha_return;

    Intrix iPeriods = Load_Intrix(Proccessed_FilePath+"iPeriods.txt", -1);
    vector<Intrix> Era_List = convertToThreeDimVec(iPeriods, n_periods, true);

    mkdir("Output_Data");
    folderName = "Output_Data/" + folderName;
    mkdir(folderName.c_str());
    for (int i = 0; i < Era_List.size(); ++i)
    {
        for (auto & iPeriod : Era_List[i])
        {
            //Calculate {beta, alpha, stock_return_akk, PERMNO, sp500_return_akk, riskFree_Return_akk}
            beta_alpha_return = Beta_Alpha_Calculate(DR_ny, iDates, sp500, riskFree, Dates, iPeriod, 130);

            //Store Data
            beta[i].insert(beta[i].end(), beta_alpha_return[0].begin(), beta_alpha_return[0].end());
            alpha[i].insert(alpha[i].end(), beta_alpha_return[1].begin(), beta_alpha_return[1].end());
            akk_return[i].insert(akk_return[i].end(), beta_alpha_return[2].begin(), beta_alpha_return[2].end());
            PERMNO[i].insert(PERMNO[i].end(), beta_alpha_return[3].begin(), beta_alpha_return[3].end());
            akk_sp500[i].insert(akk_sp500[i].end(), beta_alpha_return[4].begin(), beta_alpha_return[4].end());
            akk_riskFree[i].insert(akk_riskFree[i].end(), beta_alpha_return[5].begin(), beta_alpha_return[5].end());
        }
        // Creating a directory for the Era_List
        string dirName = folderName + "/Era";
        string periodDirName = dirName + "_" + to_string(i + 1);
        mkdir(periodDirName.c_str());

        //Save Data
        Save_Vector(periodDirName + "/beta.txt",beta[i]);
        Save_Vector(periodDirName + "/alpha.txt",alpha[i]);
        Save_Vector(periodDirName + "/akk_return.txt",akk_return[i]);
        Save_Vector(periodDirName + "/PERMNO.txt",PERMNO[i]);
        Save_Vector(periodDirName + "/akk_sp500.txt", akk_sp500[i]);
        Save_Vector(periodDirName + "/akk_riskFree.txt", akk_riskFree[i]);
    }
}
