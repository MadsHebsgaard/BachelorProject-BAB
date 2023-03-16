#include <vector>
#include <iostream>
#include <fstream>
#include "load.h"
#include <filesystem> // Requires C++17 or later //Might introduce problems for some computers
#include "calculations.h"
#include <sys/stat.h>       // For mkdir()

#include "load_output.h"



using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;
using Tensor3 = vector<Matrix>;
using Tensor4 = vector<Tensor3>;
using Tensor5 = vector<Tensor4>;


#pragma once


void BackTesting(string folderName, string dataFolderName)
{
    vector <string> infoFiles = {"/beta.txt", "/akk_return.txt", "/alpha.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/PERMNO.txt"};
    Tensor5 Data = Load_Era_Period_PrePost(dataFolderName, infoFiles);

    double betaLow = 0.5;
    double betaHigh = 1.25;

    //todo: Max_high_Beta krav?
    //todo: Min_low_Beta krav?

    int stockInf = stockInf = infoFiles.size() - 1; //Permno is excluded
    int iBeta = 0;
    int iReturn = 1;
    int iAlpha = 2;
    int iSP500 = 3;
    int iRiskFree = 4;
    int iCount = 5;
    int eraCount = Data.size();
    int maxPeriodCount = Data[0].size();

    //int periodCount = Data[eraNr].size() -1;    //TODO: Mangler 2022 sp500 data


    vector<Matrix> lowPortfolio(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));
    vector<Matrix> highPortfolio(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));
    vector<Matrix> BAB(eraCount, Matrix (maxPeriodCount, Vector(stockInf, 0)));

    vector<Matrix> preBeta(Data.size(), Matrix (maxPeriodCount, Vector(2, 0)));
    vector<Matrix> Count(Data.size(), Matrix (maxPeriodCount, Vector(2, 0)));

    Matrix BAB_Average(eraCount, Vector(stockInf+1, 0));

    //TODO: Creating the BAB portfolio as if it was a stock (ETF) with daily returns
        //TODO: plot the BAB performance VS sp500 performance

    for (int eraNr = 0; eraNr < eraCount; ++eraNr) {
        int periodCount = Data[eraNr].size();   //Todo: removed -1
        Data[eraNr].resize(periodCount);

        for (int periodNr = 0; periodNr < periodCount; ++periodNr) {
            for (int stockNr = 0; stockNr < Data[eraNr][periodNr][0][0].size(); ++stockNr) {

                double hist_Beta = Data[eraNr][periodNr][0][0][stockNr];
                //Low beta portefolio
                if(hist_Beta < betaLow)
                {
                    double ammount = betaLow - hist_Beta;
                    ammount = 1;    //ammount = 1 for equal weighted
                    preBeta[eraNr][periodNr][0] += ammount * hist_Beta;
                    Count[eraNr][periodNr][0] += ammount;
                    for (int inform = 0; inform < stockInf; ++inform) {
                        lowPortfolio[eraNr][periodNr][inform] += ammount * Data[eraNr][periodNr][1][inform][stockNr];
                    }
                }
                    //High beta portefolio
                else if(hist_Beta > betaHigh)
                {
                    double ammount = hist_Beta - betaHigh;  //ammount = 1 for equal weighted
                    ammount = 1;    //ammount = 1 for equal weighted
                    preBeta[eraNr][periodNr][1] += ammount * hist_Beta;
                    Count[eraNr][periodNr][1] += ammount;
                    for (int inform = 0; inform < stockInf; ++inform) {
                        highPortfolio[eraNr][periodNr][inform] += ammount * Data[eraNr][periodNr][1][inform][stockNr];
                    }
                }
            }
            preBeta[eraNr][periodNr][0] /= Count[eraNr][periodNr][0];
            preBeta[eraNr][periodNr][1] /= Count[eraNr][periodNr][1];

            for (int info = 0; info < stockInf; ++info) {
                lowPortfolio[eraNr][periodNr][info] /= Count[eraNr][periodNr][0];     //Average
                highPortfolio[eraNr][periodNr][info] /= Count[eraNr][periodNr][1];    //Average
            }
            //double LowAmmount = preBeta[periodNr][1] / (preBeta[periodNr][0]+preBeta[periodNr][1]);
            double BhighLow = preBeta[eraNr][periodNr][1] - preBeta[eraNr][periodNr][0];
            double lowAmmount = ( preBeta[eraNr][periodNr][1] / BhighLow);
            double highAmmount = -( preBeta[eraNr][periodNr][0] / BhighLow);

            for (int info = 0; info < stockInf; ++info) {
                BAB[eraNr][periodNr][info] = lowAmmount * lowPortfolio[eraNr][periodNr][info] + highAmmount * highPortfolio[eraNr][periodNr][info];
                BAB_Average[eraNr][info] += BAB[eraNr][periodNr][info];
            }
            BAB_Average[eraNr][iCount] += Count[eraNr][periodNr][0] + Count[eraNr][periodNr][1];
        }
        for (int info = 0; info < stockInf+1; ++info)   BAB_Average[eraNr][info] /= periodCount;
    }



    //Save BAB, low and high portfolio
    mkdir("Data/BackTest");
    string folderPath = "Data/BackTest/" + folderName;
    std::filesystem::remove_all(folderPath.c_str());
    mkdir(folderPath.c_str());
    vector<string> backTestDirs {"/BAB portfolio", "/LowBeta portfolio", "/HighBeta portfolio"};

    for (int eraNr = 0; eraNr < eraCount; ++eraNr) {
        string eraDir = folderPath + "/Era_" + to_string(eraNr+1);
        mkdir(eraDir.c_str());

        vector<string> btDirs_copy = backTestDirs;
        for(auto& dir:btDirs_copy)
        {
            dir = eraDir + dir;
            mkdir(dir.c_str());
        }
        //BAB portfolio
         for (int info = 0; info < stockInf; ++info) {
            Matrix BAB_Data(3,Vector(Data[eraNr].size()));

            for (int periodNr = 0; periodNr < Data[eraNr].size(); ++periodNr) {
                BAB_Data[0][periodNr] =           BAB[eraNr][periodNr][info];
                BAB_Data[1][periodNr] =  lowPortfolio[eraNr][periodNr][info];
                BAB_Data[2][periodNr] = highPortfolio[eraNr][periodNr][info];
            }
             for (int dirNr = 0; dirNr < btDirs_copy.size(); ++dirNr) {
                 string infoDir = btDirs_copy[dirNr] + infoFiles[info];
                 Save(infoDir, BAB_Data[dirNr]);
             }
        }
    }
    //TODO: make log file
    return;
}