#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;


void Print(Matrix A)
{
    int m, n;
    m=A.size();
    n=A[0].size();

    for(int i=0; i<m; i++)
    {
        cout << "[";
        for(int j=0; j<n;j++)
        {
            cout << setw(10) << A[i][j];
        }
        cout << "  ]" << endl;
    }
    cout << endl << endl;
}
void Print_Matrix_Transposed(Matrix A)
{
    int m, n;
    n = A.size();
    m = A[0].size();

    for(int i=0; i<m; i++)
    {
        cout << "[";
        for(int j=0; j<n; j++)
        {
            cout << setw(10) << A[j][i];
        }
        cout << "  ]" << endl;
    }
    cout << endl << endl;
}
void Matrix_overview(Matrix DR)
{
    cout << "Dim DR = " << DR.size() << endl;
    for (int i = 0; i < DR.size(); ++i)
    {
        cout << "Dim DR_" << i << " = " << DR[i].size() << endl;
    }
    //cout << DR[1][DR[1].size()-1];
}
void Print(Vector v)
{
    int m=v.size();
    for(int i=0; i<m; i++)  cout << "[" << setw(12) << v[i] << "  ]" << endl;
    cout << endl;
}
void Print(Intor v)
{
    int m=v.size();
    for(int i=0; i<m; i++)  cout << "[" << setw(12) << v[i] << "  ]" << endl;
    cout << endl;
}
