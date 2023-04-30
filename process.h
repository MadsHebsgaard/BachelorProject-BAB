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

void whatToLoad(string& Dly_Mly_Both, bool& Dly, bool& Mly)
{
    if     (Dly_Mly_Both=="Dly" || Dly_Mly_Both=="dly" )    Dly=true;
    else if(Dly_Mly_Both=="Mly" || Dly_Mly_Both=="mly" )    Mly=true;
    else if(Dly_Mly_Both=="Both"|| Dly_Mly_Both=="both")   {Dly=true; Mly=true;}
}

void Process_Files(string Dly_Mly_Both, bool Reload_Everything, bool Reload_smallFiles) //todo: DR --> Rs  *everywhere*
{
    bool Dly, Mly;  whatToLoad(Dly_Mly_Both, Dly, Mly);

    //Exo dir and required files
    string Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly;
    defineFilePaths(Dly_Mly_Both, Exo_FilePath, Proccessed_FilePath, Proccessed_FilePath_incr, Proccessed_Dly, Proccessed_Mly, Proccessed_Yly);


    vector<string> filenames = {"Dly_sp500.txt","Dly_DateList.txt", "Dly_Yly_RFR.txt", "Mly_pMC.txt", "First_MC.txt", "Yly_Inflation.txt"};
    if(Dly) filenames.push_back("Dly_Rs.txt");
    if(Mly) filenames.push_back("Mly_Rs.txt");

    if (areFilesExistInDirectory(filenames, Exo_FilePath)) {
        Intrix Rs_Dates;


        //Make dirs
        mkdir(Proccessed_FilePath.c_str());
        mkdir(Proccessed_Dly.c_str());
        mkdir(Proccessed_Mly.c_str());
        mkdir(Proccessed_Yly.c_str());

        cout << "\nCreating all necessary files and putting them in sub-directorys inside the directory " << Proccessed_FilePath << ":\n\n";

        //What to create files for:
        int max = 99999999;

        if(Dly)
        {
            //Create Daily big files:
            if (Reload_Everything || !areFilesExistInDirectory({"Rs.txt"}, Proccessed_Dly)) {
                Matrix Rs = Load_DR(Exo_FilePath+"Dly_Rs.txt", max);  //TODO: uncomment
                Compress_DR(Proccessed_Dly+"Rs.txt", Rs);   //TODO: uncomment
                cout << "Created " << Proccessed_Dly << "Rs.txt\n";
            }
            if (Reload_Everything || !areFilesExistInDirectory({"Rs_Dates.txt"}, Proccessed_Dly)) {
                Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
                Save(Proccessed_Dly+"Rs_Dates.txt", Rs_Dates);
                cout << "Created " << Proccessed_Dly << "Rs_Dates.txt\n";
            }
            else {
                Rs_Dates = Load_Intrix(Proccessed_Dly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }
            //Create Dly small files
            vector<string> smallFiles = {"Rs_iDates.txt", "iPeriods.txt", "riskFreeReturn.txt"};
            if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(smallFiles, Proccessed_Dly))
            {
                //DateList
                Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
                Dly_DateList.push_back(55555555);
                for(int i=0; i<20; i++) Dly_DateList.push_back(999999999);
                Save(Proccessed_Dly + "DateList.txt", Dly_DateList);

                //Rs_iDates
                Intrix Rs_iDates = Dly_Dates_to_iDates(Rs_Dates, Dly_DateList, 1);
                Save(Proccessed_Dly + "Rs_iDates.txt", Rs_iDates);

                //iPeriods
                Intrix iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);
                Save(Proccessed_Dly + "iPeriods.txt", iPeriods);

                //Risk free return
                Vector Dly_YearlyRFR = Load_Vector(Exo_FilePath + "Dly_Yly_RFR.txt");
                Vector Dly_RFR = DailyYearly_to_DailyDaily_Return(Dly_YearlyRFR, iPeriods);
                Save(Proccessed_Dly + "riskFreeReturn.txt", Dly_RFR);

                //sp500
                Vector sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
                Save(Proccessed_Dly+"sp500.txt",sp500);

                cout << "Created " << Proccessed_Dly << "Rs_iDates.txt\n";
                cout << "Created " << Proccessed_Dly << "iPeriods.txt\n";
                cout << "Created " << Proccessed_Dly << "riskFreeReturn.txt\n";
            }
        }
        if(Mly)
        {
            //Create Montly big files:
            if (Reload_Everything || !areFilesExistInDirectory({"Rs.txt"}, Proccessed_Mly)) {
                Matrix Rs = Load_DR(Exo_FilePath+"Mly_Rs.txt", max);  //TODO: uncomment
                Compress_DR(Proccessed_Mly+"Rs.txt", Rs);   //TODO: uncomment
                cout << "Created " << Proccessed_Mly << "Rs.txt\n";
            }
            if (Reload_Everything || !areFilesExistInDirectory({"Rs_Dates.txt"}, Proccessed_Mly)) {
                Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Mly_Rs.txt");
                Save(Proccessed_Mly+"Rs_Dates.txt", Rs_Dates);
                cout << "Created " << Proccessed_Mly << "Rs_Dates.txt\n";
            }
            else {
                Rs_Dates = Load_Intrix(Proccessed_Mly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }
            //Create Mly small files
            vector<string> smallFiles = {"Rs_iDates.txt", "iPeriods.txt", "riskFreeReturn.txt"};
            if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(smallFiles, Proccessed_Mly))
            {
                //DateList
                Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
                Intor Mly_DateList = Dly_to_Mly_DateList(Dly_DateList);
                Mly_DateList.push_back(55555555);
                for(int i=0; i<20; i++) Mly_DateList.push_back(999999999);
                Save(Proccessed_Mly + "DateList.txt", Mly_DateList);

                //Rs_iDates
                Intrix Rs_iDates = Mly_Dates_to_iDates(Rs_Dates, 1);  //TODO: Exact EOM date to iEOM function
                Save(Proccessed_Mly + "Rs_iDates.txt", Rs_iDates);

                //Mly_iPeriods
                Intrix Mly_iPeriods = x_iPeriods("Yly", "Mly", Mly_DateList);   //0-11, 12-23 , ...
                Save(Proccessed_Mly + "iPeriods.txt", Mly_iPeriods);

                //iPeriods needed to create 'Risk free return'
                Intrix Dly_M_iPeriods = x_iPeriods("Mly", "Dly", Dly_DateList);
                Intrix Dly_Y_iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);

                //Risk free return
                Vector Dly_YearlyRFR = Load_Vector(Exo_FilePath + "Dly_Yly_RFR.txt");
                Vector Dly_RFR = DailyYearly_to_DailyDaily_Return(Dly_YearlyRFR, Dly_Y_iPeriods);
                Vector Mly_RFR = period_accumulate_of_Dly(Dly_RFR, Dly_M_iPeriods);
                Save(Proccessed_Mly + "riskFreeReturn.txt", Mly_RFR);

                //sp500
                Vector sp500 = Load_Vector(Exo_FilePath+"Dly_sp500.txt");
                Vector Mly_sp500 = period_accumulate_of_Dly(sp500, Dly_M_iPeriods);
                Save(Proccessed_Mly+"sp500.txt",Mly_sp500);
            }
            cout << "Created " << Proccessed_Mly << "Rs_iDates.txt\n";
            cout << "Created " << Proccessed_Mly << "iPeriods.txt\n";
            cout << "Created " << Proccessed_Mly << "riskFreeReturn.txt\n";
        }
        //Yly files, always make these if missing or small/big reload
        vector<string> Yly_smallFiles = {"pMC.txt", "Inflation_Factor.txt"};
        if(Reload_Everything || Reload_smallFiles || !areFilesExistInDirectory(Yly_smallFiles, Proccessed_Yly))
        {
            //Rs_Dates is needed for correcting Market Cap's otherwise 2 missing stocks
            if (Dly_Mly_Both != "Dly" || Dly_Mly_Both != "dly")
            {
                if (!filePath_exists(Proccessed_Dly+"Rs_Dates.txt"))
                {
                    Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
                    Save(Proccessed_Dly+"Rs_Dates.txt", Rs_Dates);
                    cout << "Created " << Proccessed_Dly << "Rs_Dates.txt\n";
                }
                else    Rs_Dates = Load_Intrix(Proccessed_Dly+"Rs_Dates.txt", -1); //if Rs_Dates already exists
            }

            //Market Cap
            int factor = 1; //todo: needs to be 1 atm (Compress_MC)
            Matrix Mly_pMC = Load_Mly_MarketCap(Exo_FilePath+"Mly_pMC.txt", factor);
            Matrix First_MC = Load_Mly_MarketCap(Exo_FilePath + "First_MC.txt", factor);
            Matrix pMC_2comb = combine_First_with_Mly_MC(Mly_pMC, First_MC);
            Matrix pMC_3comb = Missing_MC_to_Zero(pMC_2comb, Rs_Dates);
            //Compress_MC(Proccessed_Mly + "pMC.txt", pMC_3comb, 1);
            Matrix Yly_pMC = MarketCap_Mly_to_Yly(pMC_3comb);
            Compress_MC(Proccessed_Yly + "pMC.txt", Yly_pMC, 1);


            //Inflation factors
            Vector Inflation = Load_Vector(Exo_FilePath + "Yly_Inflation.txt");
            Vector Inflation_factors = Inflation_Factors_from_yrly_inf(Inflation);
            Save(Proccessed_Yly +"Inflation_Factor.txt", Inflation_factors);
        }
    }
    else
        HowToGetStarted();
}



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

    //Dly
        //Rs
        Matrix Rs = Load_DR(Exo_FilePath+"Dly_Rs.txt", max);  //TODO: uncomment
        Compress_DR(Proccessed_Dly+"Rs.txt", Rs);   //TODO: uncomment

        //RS_Dates
        Intrix Dly_Rs_Dates = Load_Dates_from_Rs(Exo_FilePath+"Dly_Rs.txt");
        Save(Proccessed_Dly+"Rs_Dates.txt", Dly_Rs_Dates);

        //DateList
        Intor Dly_DateList = Load_Intor(Exo_FilePath + "Dly_DateList.txt");
        Dly_DateList.push_back(55555555);
        for(int i=0; i<20; i++) Dly_DateList.push_back(999999999);
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
        Mly_DateList.push_back(55555555);
        for(int i=0; i<20; i++) Mly_DateList.push_back(999999999);
        Save(Proccessed_Mly + "DateList.txt", Mly_DateList);

        //Rs_iDates
        Intrix Mly_Rs_iDates = Mly_Dates_to_iDates(Mly_Rs_Dates, 1);
        Save(Proccessed_Mly + "Rs_iDates.txt", Mly_Rs_iDates);

        //Mly_iPeriods
        Intrix Mly_iPeriods = x_iPeriods("Yly", "Mly", Mly_DateList);   //0-11, 12-23 , ...
        Save(Proccessed_Mly + "iPeriods.txt", Mly_iPeriods);

        //iPeriods needed to create 'Risk free return'
        Intrix Dly_M_iPeriods = x_iPeriods("Mly", "Dly", Dly_DateList);
        Intrix Dly_Y_iPeriods = x_iPeriods("Yly", "Dly", Dly_DateList);

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