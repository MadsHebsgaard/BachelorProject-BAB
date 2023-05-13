#include <vector>
#include <iostream>
#include <fstream>
#include "calculations.h"
#include "save.h"
#include <sys/stat.h>       // For mkdir()

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

#pragma once



void Setup_all_files()
{
    //Dir paths
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly, incr;
    defineFilePaths(incr, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);

    //Required files
    vector<string> filenames = {"Dly_sp500.txt","Dly_DateList.txt", "Dly_Yly_RFR.txt", "Mly_pMC.txt", "First_MC.txt", "Dly_Rs.txt", "Mly_Rs.txt", "Yly_Inflation.txt"};

    if (areFilesExistInDirectory(filenames, Exo_FilePath))
    {
        //Make dirs
        mkdir(Proccessed_FilePath.c_str());
        mkdir(Proccessed_Dly.c_str());
        mkdir(Proccessed_Mly.c_str());
        mkdir(Proccessed_Yly.c_str());

        cout << "\nCreating all necessary files and putting them in sub-directorys inside the directory " << Proccessed_FilePath << ":\n";
        cout << "This will likely take 10 to 25 minutes, as large files are being read and new files are being calculated and compressed.\n\n";
        cout << "The function only needs to be activated once.\n\n";
        int max = 99999999;

        //When testing to reduce load time by 99%
        //Intrix Dly_Rs_Dates = Load_Intrix("Data/Input/Processed_Files/Dly/Rs_Dates.txt",-1);
        //Save(Proccessed_Dly+"Rs_Dates.txt", Dly_Rs_Dates);
        //Intrix Mly_Rs_Dates = Load_Intrix("Data/Input/Processed_Files/Mly/Rs_Dates.txt",-1);
        //Save(Proccessed_Mly+"Rs_Dates.txt", Mly_Rs_Dates);

    //Dly
        //Rs
        Matrix Rs = Load_DR(Exo_FilePath+"Dly_Rs.txt", max);  //TODO: uncomment
        Compress_DR(Proccessed_Dly+"Rs.txt", Rs);   //TODO: uncomment

        //RS_Dates
        Intrix Dly_Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
        Save(Proccessed_Dly+"Rs_Dates.txt", Dly_Rs_Dates);

        //DateList
        Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
        Save(Proccessed_Dly + "DateList.txt", Dly_DateList);

        //Rs_iDates
        Intrix Dly_Rs_iDates = Dly_Dates_to_iDates(Dly_Rs_Dates, Dly_DateList, 1);
        Save(Proccessed_Dly + "Rs_iDates.txt", Dly_Rs_iDates);

        //iPeriods
        Intrix Dly_iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);
        Save(Proccessed_Dly + "iPeriods.txt", Dly_iPeriods);

        //Risk free return
        Vector Dly_YearlyRFR = Load_Vector(Exo_FilePath + "Dly_Yly_RFR.txt");
        Vector Dly_RFR = DailyYearly_to_DailyDaily_Return(Dly_YearlyRFR, Dly_iPeriods);
        Save(Proccessed_Dly + "riskFreeReturn.txt", Dly_RFR);

        //sp500
        Vector Dly_sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
        Save(Proccessed_Dly+"sp500.txt", Dly_sp500);


    //Mly
        //Rs
        Matrix Mly_Rs = Load_DR(Exo_FilePath+"Mly_Rs.txt", max);
        Compress_DR(Proccessed_Mly+"Rs.txt", Mly_Rs);

        //RS_Dates
        Intrix Mly_Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Mly_Rs.txt");
        Save(Proccessed_Mly+"Rs_Dates.txt", Mly_Rs_Dates);

        //DateList
        Intor Mly_DateList = Dly_to_Mly_DateList(Dly_DateList);
        Save(Proccessed_Mly + "DateList.txt", Mly_DateList);

        //Rs_iDates
        Intrix Mly_Rs_iDates = Mly_Dates_to_iDates(Mly_Rs_Dates, 1);
        Save(Proccessed_Mly + "Rs_iDates.txt", Mly_Rs_iDates);

        //Mly_iPeriods
        //Intrix Mly_iPeriods = x_iPeriods("Yly", "Mly", Mly_DateList);   //0-11, 12-23 , ...
        Intrix Mly_iPeriods = Create_Mly_iPeriods(Mly_DateList);   //0-11, 12-23 , ...
        Save(Proccessed_Mly + "iPeriods.txt", Mly_iPeriods);

        //iPeriods needed to create 'Risk free return'
        Intrix Dly_M_iPeriods = x_iPeriods("Mly", "Dly", Dly_DateList);

        //Risk free return
        Vector Mly_RFR = period_accumulate_of_Dly(Dly_RFR, Dly_M_iPeriods);
        Save(Proccessed_Mly + "riskFreeReturn.txt", Mly_RFR);

        //sp500
        Vector Mly_sp500 = period_accumulate_of_Dly(Dly_sp500, Dly_M_iPeriods);
        Save(Proccessed_Mly+"sp500.txt", Mly_sp500);


    //Yly
        //Market Cap
        int factor = 1; //todo: needs to be 1 atm (Compress_MC)
        Matrix Mly_pMC = Load_Mly_MarketCap(Exo_FilePath+"Mly_pMC.txt", factor);
        Matrix First_MC = Load_Mly_MarketCap(Exo_FilePath + "First_MC.txt", factor);
        Matrix pMC_2comb = combine_First_with_Mly_MC(Mly_pMC, First_MC);
        Matrix pMC_3comb = Missing_MC_to_Zero(pMC_2comb, Dly_Rs_Dates);
        //Compress_MC(Proccessed_Mly + "pMC.txt", pMC_3comb, 1);    //If monthly was of use
        Matrix Yly_pMC = MarketCap_Mly_to_Yly(pMC_3comb);
        Compress_MC(Proccessed_Yly + "pMC.txt", Yly_pMC, 1);

        //Inflation factors
        Vector Inflation = Load_Vector(Exo_FilePath + "Yly_Inflation.txt");
        Vector Inflation_factors = Inflation_Factors_from_yrly_inf(Inflation);
        Save(Proccessed_Yly +"Inflation_Factor.txt", Inflation_factors);
    }
    else
        HowToGetStarted();
}