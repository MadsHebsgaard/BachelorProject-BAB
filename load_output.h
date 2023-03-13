

#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "calculations.h"

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;


vector<Matrix> Load_Era_PrePost(string runName, int Era_nr)
{
    string dirPath = "Data/Output/Double/" + runName + "/Era_"+to_string(Era_nr);
    vector<string> periodFolders = {dirPath+"/Pre_Period", dirPath+"/Period"};
    vector<string> fileNames = {"/PERMNO.txt", "/beta.txt", "/alpha.txt", "/akk_return.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
    vector<Matrix> Data(2, Matrix(fileNames.size(),Vector(0)));
    double number = 1;

    for (int prePost = 0; prePost < 2; ++prePost) {
        for (int fileNr = 0; fileNr < fileNames.size(); ++fileNr) {
            ifstream fil(periodFolders[prePost]+fileNames[fileNr]);
            while(fil >> number) {
                Data[prePost][fileNr].push_back(number);
            }
        }
    }
    return Data;
}
vector<vector<Matrix>> Load_Era_n_PrePost_Period(string runName, int Era_nr)
{
    string dirPath = "Data/Output/Era_PrePost_Period/" + runName + "/Era_"+to_string(Era_nr);
    vector<string> periodFolders = {dirPath+"/Pre_Period", dirPath+"/Period"};
    vector<string> fileNames = {"/PERMNO.txt", "/beta.txt", "/alpha.txt", "/akk_return.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};

    int Period_count = numSubdirsInDir(periodFolders[0]);
    vector<vector<Matrix>> Data(2, vector<Matrix>(Period_count, Matrix(fileNames.size(),Vector(0))));
    double number;

    for (int prePost = 0; prePost < 2; ++prePost) {
        for (int Period = 0; Period < Period_count; ++Period) {
            string periodPath = periodFolders[prePost] + "/Period_" + to_string(Period+1);
            for (int fileNr = 0; fileNr < fileNames.size(); ++fileNr) {
                ifstream fil(periodPath + fileNames[fileNr]);
                while (fil >> number) {
                    Data[prePost][Period][fileNr].push_back(number);
                }
            }
        }
    }
    return Data;
}
vector<vector<vector<Matrix>>> Load_Era_PrePost_Period(string runName)
{
    //vector<string> fileNames = {"/beta.txt", "/akk_return.txt", "/alpha.txt", "/PERMNO.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
    vector<string> fileNames = {"/beta.txt", "/akk_return.txt", "/alpha.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/PERMNO.txt"};


    vector<string> prePostName = {"/Pre_Period", "/Period"};
    string eraPaths = "Data/Output/Era_PrePost_Period/" + runName;
    int Era_count = numSubdirsInDir(eraPaths);
    double number;
    int maxEraPeriods = numSubdirsInDir(eraPaths+"/Era_1/Period");

    vector<vector<vector<Matrix>>> Data(Era_count, vector<vector<Matrix>>(2, vector<Matrix>(maxEraPeriods, Matrix(fileNames.size(),Vector(0)))));

    for (int eraNr = 0; eraNr < Era_count; ++eraNr) {
        string thisEraPath = eraPaths + "/Era_" + to_string(eraNr+1);

        for (int prePost = 0; prePost < 2; ++prePost) {
            string prePostPath = thisEraPath + prePostName[prePost];
            int Period_count = numSubdirsInDir(prePostPath);

            Data[eraNr][prePost].resize(Period_count);

            for (int Period = 0; Period < Period_count; ++Period) {
                string periodPath = prePostPath + "/Period_" + to_string(Period+1);

                for (int fileNr = 0; fileNr < fileNames.size(); ++fileNr) {
                    ifstream fil(periodPath + fileNames[fileNr]);
                    while (fil >> number) {
                        Data[eraNr][prePost][Period][fileNr].push_back(number);
                    }
                }
            }
        }
    }
    return Data;
}
vector<vector<vector<Matrix>>> Load_Era_Period_PrePost(string runName)
{

    //vector<string> fileNames = {"/PERMNO.txt", "/beta.txt", "/alpha.txt", "/akk_return.txt", "/akk_sp500.txt", "/akk_riskFree.txt"};
    vector<string> fileNames = {"/beta.txt", "/akk_return.txt", "/alpha.txt", "/akk_sp500.txt", "/akk_riskFree.txt", "/PERMNO.txt"};

    vector<string> prePostName = {"/Pre_Period", "/Period"};
    string eraPaths = "Data/Output/Era_Period_PrePost/" + runName;
    int Era_count = numSubdirsInDir(eraPaths);
    double number;
    int maxEraPeriods = numSubdirsInDir(eraPaths+"/Era_1");

    vector<vector<vector<Matrix>>> Data(Era_count, vector<vector<Matrix>>(maxEraPeriods, vector<Matrix>(2, Matrix(fileNames.size(),Vector(0)))));

    for (int eraNr = 0; eraNr < Era_count; ++eraNr) {
        string thisEraPath = eraPaths + "/Era_" + to_string(eraNr+1);
        int Period_count = numSubdirsInDir(thisEraPath);
        Data[eraNr].resize(Period_count);

        for (int Period = 0; Period < Period_count; ++Period) {
            string periodPath = thisEraPath + "/Period_" + to_string(Period+1);

        for (int prePost = 0; prePost < 2; ++prePost) {
            string prePostPath = periodPath + prePostName[prePost];

                for (int fileNr = 0; fileNr < fileNames.size(); ++fileNr) {
                    ifstream fil(prePostPath + fileNames[fileNr]);
                    while (fil >> number) {
                        Data[eraNr][Period][prePost][fileNr].push_back(number);
                    }
                }
            }
        }
    }
    return Data;
}
