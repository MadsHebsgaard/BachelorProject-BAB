#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;
using Intor = vector<int>;
using Intrix = vector<vector<int>>;
using Vector = vector<double>;
using Matrix = vector<vector<double>>;

    //Essential
Matrix Alpha_Beta_Mu(Matrix DR, Intrix iDates, Vector sp500, const Intor& iPeriod, int minTradingDaysInPeriod); //TODO
double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR);      //TODO
double Calculate_stockReturn(Vector stock, int iLength, int iStart_DR);             //TODO
double CovSP500(Vector Stock, Vector sp500, int iStart, int iEnd, int iStart_DR);   //TODO: Fix start date and end date
double Var(Vector sp500, int Start_number, int End_number);                         //TODO: Fix start date and end date

    //temporary
void day_est_checker();
int day_after1926(int date);
int Find_date_integer(int date, Intor& Date_list);
void day_est_checker(int max);
void Date_Checker(int max);
double Vector_Average(const Vector& v);

    //Support
Vector Column_of_Matrix(Matrix A, int n);
Intor Column_of_Intrix(Intrix A, int n);
Intrix Dates_to_iDates(Intrix Dates, Intor True_Dates, int start_j);

    //UI
void Print_Matrix(Matrix A);
void Print_Matrix_Transposed(Matrix A);
void Print_Vector(Vector v);
void Print_Intor(Intor v);
void Matrix_overview(Matrix DR);

    //Load from fil
Matrix Load_Matrix(const string& fn);
Matrix Load_DR(const string& fn, int max);
Intrix Load_Dates_from_DR(const string& fn);
Matrix Load_DR_Compressed(const string& fn, int max);
Intrix Load_StockDays_Compressed(const string& fn, bool With_Id, int max);
Intrix Load_StockDays_from_DR(const string& fn, int max);
Intrix Load_Intrix(const string& fn, int max);
Intrix Load_Dates(const string& fn);
Intor Load_Intor(const string& fn);
Vector Load_Vector(const string& fn);

    //Save to file
void Save_Vector(const string& fn, const Vector& v);
void Save_Intor(const string& fn, const Intor& v);
void Save_Intrix(const string& fn, Intrix A);
void Compress_DR(const string& fn, Matrix DR);
void Compress_DR_StockDays(const string& fn, Intrix StockDays);
void sp500_Dates_to_monthly_dDate_periods();

Matrix Edit_DR(Matrix A, double max_ratio, int minTradingDays);
Vector Matrix_Column(Matrix A, int j);
Intrix Remove_Missing_ID_Intrix(Intrix A, Vector v);

//Intrix sp500_Dates_to_monthly_dDate_periods(Intor allDates)

void Create_Files_from_DR(int max);
void TestCalculations(int max);
void RunCalculations(int max);

int main()
{
    //sp500_Dates_to_monthly_dDate_periods();
    //Create_Files_from_DR(-1);
    TestCalculations(5000);
    //RunCalculations(500);
    return 0;
}

void TestCalculations(int max)
{
    //max = 500;
    int minTradingDays = 75;

        //DR with condition for inclusion
    //Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    //Matrix DR_ny = Edit_DR(DR,0.3,minTradingDays);

    Matrix DR_ny = Load_DR_Compressed("DR_Compressed.txt", max);
    //Intrix StockDays = Load_StockDays_Compressed("DR_StockDays.txt", true, max);
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    //Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);
    Intrix iPeriod = Load_Intrix("Yearly_iPeriods.txt", -1);
    Vector sp500 = Load_Vector("sp500.txt");
    Intor Dates = Load_Intor("SP_Dates.txt");


    //iDates with same stocks as DR_ny
    //Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    //iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

        //Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    //Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);

    Vector beta(0); //TODO: One for mu and alpha too

    //Beta
    for (int i = 0; i < iPeriod.size(); ++i)
    {
        Matrix beta_mean = Alpha_Beta_Mu(DR_ny, iDates, sp500, iPeriod[i], 130);
        Data.push_back(beta_mean);

        if(i < 15)
        {
            cout << setw(9) << "n = " << beta_mean[0].size();
            cout << setw(18) << "Avg_B = " << Vector_Average(beta_mean[0]) << endl;
        }
        beta.insert(beta.end(), beta_mean[0].begin(), beta_mean[0].end());
    }
    //cout << "mean beta = " << Vector_Average(Data[0]) << endl;
}
void RunCalculations(int max)
{
    //max = 5000;
    int minTradingDays = 75;

    //DR with condition for inclusion
    Matrix DR = Load_DR_Compressed("DR_Compressed.txt", max);
    Matrix DR_ny = Edit_DR(DR,0.3,minTradingDays);

    //iDates with same stocks as DR_ny
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);
    iDates = Remove_Missing_ID_Intrix(iDates, Matrix_Column(DR_ny, 0));

    //Intrix iPeriod = Load_Intrix("iPeriod.txt", -1);
    Intrix iPeriod = Load_Intrix("Monthly_iPeriods.txt", -1);

    //Load sp500 stock returns
    Vector sp500 = Load_Vector("sp500.txt");

    cout << "Files was Loaded.\n\n";

    vector<Matrix> Data(0);

    //Beta
    for (int i = 0; i < iPeriod.size(); ++i)
    {
        Matrix beta_mean = Alpha_Beta_Mu(DR_ny, iDates, sp500, iPeriod[i], 8);
        Data.push_back(beta_mean);

        cout << setw(9) << "n = " << beta_mean[0].size();
        cout << setw(18) << "Avg_B = " << Vector_Average(beta_mean[0]) << endl;
    }
}

void sp500_Dates_to_monthly_dDate_periods()
{
    Intor allDates = Load_Intor("SP_Dates.txt");
    int i = 0, periodStart = allDates[0];
    Intrix MonthlyPeriods(0);

    while (allDates[i] < 99999999)
    {
        if( (allDates[i] - periodStart) > 40)
        {
            MonthlyPeriods.push_back({periodStart, allDates[i-1]});
            periodStart = allDates[i];
        }
        i++;
    }
    //Save_Intrix("Monthly_dPeriods.txt", MonthlyPeriods);
    Intrix Monthly_iPeriods = Dates_to_iDates(MonthlyPeriods,allDates,0);
    Save_Intrix("Monthly_iPeriods.txt", Monthly_iPeriods);


    i = 0;
    periodStart = allDates[0];
    Intrix YearlyPeriods(0);
    while (allDates[i] < 99999999)
    {
        if( (allDates[i] - periodStart) > 4000)
        {
            YearlyPeriods.push_back({periodStart, allDates[i-1]});
            periodStart = allDates[i];
        }
        i++;
    }
    Intrix Yearly_iPeriods = Dates_to_iDates(YearlyPeriods,allDates,0);
    Save_Intrix("Yearly_iPeriods.txt", Yearly_iPeriods);
}
Intrix Remove_Missing_ID_Intrix(Intrix A, Vector v)
{
    Intrix B(0);
    size_t k=0;
    for (auto & stock : A)
    {
        //if(k==v.size()-1) break;
        if(v[k] == stock[0]) {k++;   B.push_back(stock);}
    }
    return B;
}
Vector Matrix_Column(Matrix A, int j)
{
    Vector col(A.size());
    for (int i = 0; i < A.size(); ++i)  col[i] = A[i][j];
    return col;
}
Matrix Edit_DR(Matrix A, double max_ratio, int minTradingDays)
{
    Matrix DR_ny(0);
    int træk, zero_total;
    int max_træk = 6;
    int replace = -2;
    int continue_amount = 0;

    for (auto & stock : A)
    {
        zero_total = 0;
        træk = 0;
        if (stock.size() < minTradingDays)   {continue_amount++; continue;}
        for (int j = 0; j < stock.size(); ++j)
        {
            if(stock[j] == 0 or stock[j] == replace)    {træk++; zero_total++;}
            else                 træk = 0;

            if      (træk == max_træk)  for(int k = 1; k <= max_træk; ++k) stock[j - max_træk + k] = replace;
            else if (træk > max_træk) stock[j] = replace;
        }
        if((zero_total+0.0) / stock.size() < max_ratio)  {DR_ny.push_back(stock);           /*cout << endl << "zero_total:" << zero_total << "  size: " << A[stock].size();*/}
        //else cout << endl << stock << "   " << A[stock][0] << "  " << (zero_total+0.0)/A[stock].size();
    }

    cout << "\ncontinue_amount = " << continue_amount;
    cout << "\nDR_ny size = " << DR_ny.size();
    cout << "\nZero ratio too high = " << A.size() - DR_ny.size() - continue_amount;
    cout << "\n\n";

    return DR_ny;
}
double Vector_Average(const Vector& v)
{
    double sum=0;
    for (double e : v)  sum+=e;
    return sum/v.size();
}
void Create_Files_from_DR(int max)
{
    //Need files:
        //DR_NO_Ticker.txt
        //SP_Dates
        //Period_Dates.txt  (optional if needing periods)

    //DR
    Matrix DR = Load_DR("DR_No_Ticker.txt", max);
    Compress_DR("DR_Compressed.txt", DR);

    //StockDays - Usefull while testing
    Intrix StockDays = Load_StockDays_from_DR("DR_No_Ticker.txt", max);
    Compress_DR_StockDays("DR_StockDays.txt", StockDays);

    //DR_Dates
    Intrix DR_Dates = Load_Dates_from_DR("DR_No_Ticker.txt");
    Save_Intrix("DR_Dates.txt", DR_Dates);

    //iDates
    Intor SP_Dates = Load_Intor("SP_Dates.txt");
    Intrix iDates = Dates_to_iDates(DR_Dates, SP_Dates, 1);
    Save_Intrix("DR_iDates.txt", iDates);

    //iPeriod
    Intrix dPeriod = Load_Intrix("Period_Dates.txt", -1);
    Save_Intrix("dPeriod.txt", dPeriod);

    //Intrix dPeriod = Load_Intrix("dPeriod.txt", -1);
    Intrix iPeriod = Dates_to_iDates(dPeriod, SP_Dates, 0);
    Save_Intrix("iPeriod.txt", iPeriod);

    sp500_Dates_to_monthly_dDate_periods();
}
void Date_Checker(int max)
{
    if (max < 0)    max = 36146;
    Intrix StockDays = Load_StockDays_Compressed("DR_StockDays.txt",true, max);
    Intor Dates = Load_Intor("SP_Dates.txt");
    Intrix iDates = Load_Intrix("DR_iDates.txt", max);

    //78094 & 15451

    cout << endl;

    for (int i = 0; i < max; ++i)
    {
        //int x = 300;   //1000  Last is unequal 19961105 19961106 78094   20070625
        for(int x = 0; x < iDates[i][2]-iDates[i][1]; ++x)
        {
            if(iDates[i][2] - iDates[i][1] >= x)
            {
                if(Dates[iDates[i][1] + x] != StockDays[i][1 + x])
                {
                    cout << "Last is unequal " << Dates[iDates[i][1] + x];
                    cout << " " << StockDays[i][1 + x] << " " << StockDays[i][0] << "   ";
                    cout << Dates[iDates[i][2]] << endl;
                }
            }
            else
            {
                //cout << setw(5) << iDates[i][2] - iDates[i][1] + 1 << " < " << x+1 << ". ID = " << iDates[i][0] << endl;
            }
        }
    }
}
double Calculate_stockBeta(Vector stock, const Vector& sp500, int iStart, int iEnd, int iStart_DR)
{
    double stockCov = CovSP500(stock, sp500, iStart, iEnd, iStart_DR);
    double stockSP500var = Var(sp500, iStart, iEnd);
    return (stockCov/stockSP500var);
}
double Calculate_stockReturn(Vector stock, int iLength, int iStart_DR)
{
    double Return=1;
    int trading_Days = 0;
    for (int i = iStart_DR; i < iStart_DR + iLength; ++i)
    {
        if(stock[i] == -2) continue;
        trading_Days++;
        Return*=(1 + stock[i]);
    }
    return Return - 1;                                    //TODO: Return type
    //return pow(Return,263.5104/trading_Days)-1;        //
}
Matrix Alpha_Beta_Mu(Matrix DR, Intrix iDates, Vector sp500, const Intor& iPeriod, int minTradingDaysInPeriod)
{
    Vector cov(0);
    Vector sp500_Var(0);
    Vector beta(0);
    Vector returns(0);
    Vector alpha(0);
    int iStart, iEnd, iStart_DR, iLength;
    double stockBeta, stockReturn, stockAlpha;
    int emptyCount;
    //double stockCov, stockSP500var;



    for (int i = 0; i < DR.size(); ++i)
    {
        //Determen index start and end in sp500 and index to start in DR[i]
        if(iDates[i][1] < iPeriod[0])
        {
            iStart = iPeriod[0];                            //sp500 date to start
            iStart_DR = 1 + (iPeriod[0] - iDates[i][1]);    //index where period start in DR[i]
        }
        else
        {
            iStart = iDates[i][1];
            iStart_DR = 1;
        }
        iEnd = min(iDates[i][2], iPeriod[1]);               //sp500 date to end
        iLength = iEnd - iStart + 1;                        //iDates total

        //Check if stock is relevant
        emptyCount=0;
        if(iLength > 0) for (int j = iStart_DR; j < iStart_DR+iLength; ++j)     if(DR[i][j] == -2)    emptyCount++;
        if((iLength - emptyCount) < minTradingDaysInPeriod) continue;

        //Calculate beta
        stockBeta = Calculate_stockBeta(DR[i], sp500, iStart, iEnd, iStart_DR);
        beta.push_back(stockBeta);

        //Calculate return
        stockReturn = Calculate_stockReturn(DR[i], iLength, iStart_DR);
        returns.push_back(stockReturn);

        double rf_return = 0.02;    //TODO: insert periods risk free return as an argument in the function Alpha_Beta_Mu
        stockAlpha = (stockReturn-rf_return)/stockBeta;
        alpha.push_back(stockAlpha);

        //cout << "stockReturn = " << setw(10) << stockReturn << " i = " << i << endl;
    }


    double sp500_mean = Calculate_stockReturn(sp500, (iPeriod[1] - iPeriod[0]) + 1, iPeriod[0]);
    //cout << setw(18) << "sp_r = " << sp500_mean;
    //cout << "\nreturns = \n";    Print_Vector(returns);
    //cout << setw(18) << "avg_r = " << Vector_Average(returns);
    //if(DR.size() < 50) {cout << "cov = \n"; Print_Vector(cov); cout << "sp500_Var = \n"; Print_Vector(sp500_Var);}
    return {beta, returns, alpha};
}
Intrix Load_StockDays_Compressed(const string& fn, bool With_Id, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }
    cout << "Loading " << max << " StockDays from compressed file\n";
    Intor days;
    Intrix StockDays;
    int n, i=0;

    while(fil >> n)
    {
        if(!With_Id)     n--;
        days = Intor(n);
        int j = 0;
        if(With_Id) {fil >> days[0];    j++;}
        while (j < n)
        {
            fil >> days[j];
            j++;
        }
        StockDays.push_back(days);
        if(i % 2000 == 0)    cout << "Loaded up to ID: " << days[0] << ", i = " << i << endl;
        i++;
        if (i == max) break;
    }
    cout << "Sucssesfully loaded " << i << " stocks to StockDays.\n";
    return StockDays;
}
void Compress_DR_StockDays(const string& fn, Intrix StockDays)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    for (auto & stock : StockDays)
    {
        fil << endl << endl << stock.size() << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock.size(); ++j)  fil << " " << stock[j];
    }
    cout << "Compress_DR_StockDays: Sucssesfully Compressed/Saved daily 'dates' of " << StockDays.size() << " stocks to " << fn << ".\n";
}
Intrix Load_StockDays_from_DR(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }
    Intrix Stockdays;
    Intor days;
    int IDnew, date, date_old, i=0;
    double r_;
    string junk;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;
    cout << "Loaded 0/36148 = 0%.\n";

    while(!fil.eof())
    {
        days = Intor(0);
        days.push_back(IDnew);

        fil >> date;
        days.push_back(date);
        fil >> r_;
        if(r_<10000)  fil >> IDnew;
        else          IDnew = r_;

        while(true)
        {
            if(IDnew != days[0]) break;
            else if(fil.eof()) break;
            date_old = date;
            fil >> date;

            while(date_old == date)
            {
                fil >> r_;
                if(r_ < 10000)  fil >> IDnew;
                else    IDnew = r_;
                if(IDnew != days[0])    break;
                date_old = date;
                fil >> date;
                if(fil.eof()) return Stockdays;
            }
            if(IDnew != days[0])    break;
            days.push_back(date);
            fil >> r_;
            if(r_ < 10000)
            {
                fil >> IDnew;
            }
            else
            {
                IDnew = r_;
            }
        }
        Stockdays.push_back(days);
        i++;
        if(i % 500 == 0) cout << "Loaded " << i << "/36148 = " << i/361.48 << "%.\n";
        if(i == max)   break;
    }
    cout << "Loaded " << i << "/36148 = " << i / 361.48 << "%.\n";
    return Stockdays;
}
Vector Load_Vector(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Vector(0);   }
    Vector v;
    double e;
    while(!fil.eof())
    {
        fil >> e;
        v.push_back(e);
    }
    return v;
}
Intrix Dates_to_iDates(Intrix Dates, Intor True_Dates, int start_j)
{
    for (auto & Date : Dates)
    {
        for (int j = start_j; j < Date.size(); ++j)
        {
            Date[j] = Find_date_integer(Date[j], True_Dates);
        }
    }
    return Dates;
}
int Find_date_integer(int date, Intor& Date_list)
{
    int index_apprx = day_after1926(date);
    while(Date_list[index_apprx] < date)    index_apprx += 15;
    while(Date_list[index_apprx] > date)    index_apprx -= 1;
    return index_apprx;
}
Intor Load_Intor(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_Intor: Could not read the file " << fn << ".";  return Intor(0);   }
    Intor v;
    int e;
    while(!fil.eof())
    {
        fil >> e;
        v.push_back(e);
    }
    cout << "Load_Intor: Sucssesfully loaded " << v.size() << " elements from " << fn << " to Intor.\n";
    return v;
}
Intrix Load_Intrix(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_Intrix: Could not read the file " << fn << ".";  return Intrix(0);   }

    if(max < 1) max = 999999;

    int m, n;
    fil >> m;
    fil >> n;
    Intrix A(min(m, max), Intor(n));

    for (int i = 0; i < min(m, max); ++i)
        for (int j = 0; j < n; ++j)
            fil >> A[i][j];

    if ( A[0].size() == A[1].size() and A[0].size() == A[A.size()-1].size() )
            cout << "Load_Intrix: Sucessfully loaded (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
    else    cout << "Load_Intrix: Sucssesfully loaded (" << A.size() << " x X) to Intrix.\n";

    return A;
}
void Save_Intrix(const string& fn, Intrix A)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    fil << A.size() << endl;
    fil << A[0].size() << endl << endl;

    for (int i = 0; i < A.size(); ++i)
    {
        //fil << setw(10) << A[i][0];
        for (int j = 0; j < A[0].size(); ++j)
        {
            fil << setw(10) << A[i][j];
        }
        fil << endl;
    }
    cout << "Save_Intrix: Sucessfully saved (" << A.size() << " x " << A[0].size() << ") Intrix to " << fn << ".\n";
}
Matrix Load_DR_Compressed(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Load_DR_Compressed: Could not read the file " << fn << ".";  return Matrix(0);   }

    Vector Stock;
    Matrix DR;
    int n, i=0;

    while(fil >> n)
    {
        Stock = Vector(n);
        for (int j = 0; j < n; ++j)
        {
            fil >> Stock[j];
        }
            DR.push_back(Stock);
            if(i % 500 == 0)    cout << "ID = " << Stock[0] << " i = " << i << endl;
            i++;
            if (i == max) break;
    }
    cout << "Load_DR_Compressed: Sucssesfully loaded " << i << " stocks to DR from " << fn << ".\n";
    return DR;
}// vector || pushback
void Compress_DR(const string& fn, Matrix DR)
{
    ofstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return;   }

    for (auto & stock : DR)
    {
        fil << endl << endl << stock.size() << " " << stock[0] << endl << stock[1];
        for (int j = 2; j < stock.size(); ++j)  fil << " " << stock[j];
    }
    cout << "Compress_DR: Sucssesfully Compressed/Saved dily 'returns' of " << DR.size() << "  stocks to " << fn << ".\n";
}
Intrix Load_Dates_from_DR(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }
    Intrix Dates;
    Intor Stock;
    int IDnew, date, i=0;
    double r;
    string junk;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;
    cout << "Loaded 0/36148 = 0%.\n";

    while(!fil.eof())
    {
        Stock = Intor(0);
        fil >> date;
        Stock.push_back(IDnew);
        Stock.push_back(date);
        while(!fil.eof())
        {
            fil >> r;
            if(r < 10000)   fil >> IDnew;
            else            IDnew = r;
            if(IDnew != Stock[0])   break;
            fil >> date;
        }
        Stock.push_back(date);
        Dates.push_back(Stock);
        i++;
        if(i % 500 == 0) cout << "Loaded " << i << "/36148 = " << i/361.48 << "%.\n";
    }
    cout << "Loaded " << i << "/36148 = " << i / 361.48 << "%.\n";
    return Dates;
}
Matrix Load_DR(const string& fn, int max)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Matrix(0);   }
    Matrix DR;
    Vector Stock;
    int IDnew, date, date_old, i=0;
    double r, r_;
    string junk;
    int Blanc_r = -2;

    fil >> junk;
    fil >> junk;
    fil >> junk;
    fil >> IDnew;

    cout << "Loaded 0/36148 = 0%.\n";

    while(!fil.eof())
    {
        Stock = Vector(0);
        Stock.push_back(IDnew);

        {
            if(fil.eof()) break;
            fil >> date;
            fil >> r;
            if(r < 10000)   fil >> IDnew;
            else
            {
                IDnew = r;
                r = Blanc_r;     //When blanc/no return
            }
            Stock.push_back(r);
        }
        while(true)
        {
            if(IDnew != Stock[0]) break;
            else if(fil.eof()) break;
            date_old = date;
            fil >> date;

            while(date_old == date)
            {
                fil >> r;
                if(r < 10000)  fil >> IDnew;
                else    {IDnew = r; r = Blanc_r;}
                if(IDnew != Stock[0])    break;
                date_old = date;
                fil >> date;
                if(fil.eof()) break;
            }
            if(IDnew != Stock[0] or fil.eof())    break;
            fil >> r_;
            if(r_ < 10000)
            {
                r = r_;
                fil >> IDnew;
            }
            else
            {
                r = Blanc_r;
                IDnew = r_;
            }
            Stock.push_back(r);
        }
        DR.push_back(Stock);
        i++;
        if(i % 500 == 0) cout << "Loaded " << i << "/36148 = " << i/361.48 << "%.\n";
        if(i==max)   break;
    }
    cout << "Loaded " << i << "/36148 = " << i / 361.48 << "%.\n";
    return DR;
}
Intrix Load_Dates(const string& fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Intrix(0);   }

    Intrix A;
    string junk_string;
    int IDnew, junk_number, i=0;
    vector<int> Stock(3);
    fil >> junk_string;    fil >> junk_string;    fil >> junk_string;

    while(fil >> IDnew)
    {
        Stock[0] = IDnew;
        fil >> Stock[1];
        fil >> Stock[2];

        while(IDnew == Stock[0] && !fil.eof())
        {
            fil >> IDnew;
            fil >> junk_number;
            fil >> junk_number;
        }
        //cout << "ID = " << Stock[0] << ", i = " << i << endl;
        //if(i % 500 == 0) cout << "Stock[0] = " << Stock[0] << ", i = " << i << endl;
        A.push_back(Stock);
        i++;
    }
    cout << "All " << i << " dates was loaded sucsessfully.";
    return A;
}
void Print_Vector(Vector v)
{
    int m=v.size();
    for(int i=0; i<m; i++)  cout << "[" << setw(12) << v[i] << "  ]" << endl;
    cout << endl;
}
void Print_Intor(Intor v)
{
    int m=v.size();
    for(int i=0; i<m; i++)  cout << "[" << setw(12) << v[i] << "  ]" << endl;
    cout << endl;
}
double CovSP500(Vector Stock, Vector sp500, int iStart, int iEnd, int iStart_DR)
{
    int d_iDates = iEnd - iStart + 1;
    double Mean_Stock = 0, Mean_sp500 = 0, Covv = 0;
    int trading_Days = 0;

    //Mean
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i] != -2)
        {
            trading_Days++;
            Mean_Stock += Stock[i + iStart_DR];
            Mean_sp500 += sp500[i + iStart];
        }
    }
    Mean_Stock /= trading_Days;
    Mean_sp500 /= trading_Days;

        //CovSP500
    for (int i = 0; i < d_iDates; i++)
    {
        if(Stock[i] != -2)
            Covv += (Stock[i+iStart_DR] - Mean_Stock) * (sp500[i + iStart] - Mean_sp500);
    }
    return Covv /(trading_Days - 1);
}
double Var(Vector sp500, int Start_number, int End_number)
{
    double Mean = 0, Variance = 0;
    //int d_dates = End_number-Start_number+1;
    //if(sp500.size() < d_dates) d_dates = sp500.size();

    int trading_Days=0;

    //Mean
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
        {
            trading_Days++;
            Mean += sp500[i];
        }
    }
    Mean /= trading_Days;

    //Variance
    for (int i = Start_number; i < End_number+1; i++)
    {
        if(sp500[i] != -2)
            Variance += (sp500[i] - Mean) * (sp500[i] - Mean);
    }
    return Variance /(trading_Days - 1);
}
int day_after1926(int date)
{
    int year_after = date/10000 - 1926;
    int month_after = (date%10000)/100 - 1;
    int days_after = date%100-2;
    return 263 * year_after + month_after * 23 + days_after * 13/16;
} //TODO: Check if it works decently for years after 1927
void Print_Matrix(Matrix A)
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
Matrix Load_Matrix(const string& fn)
{
    //Forbind til fil
    ifstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return Matrix(0);   }

    //Lav matrix struktur
    size_t m,n;
    fil >> m;
    fil >> n;
    Matrix A(m);
    for (auto& r: A)    r.resize(n);

    //Indhent matrix elementer
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)  fil >> A[i][j];
    }

    return A;
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
Vector Column_of_Matrix(Matrix A, int n)
{
    Vector v(A.size());
    for (int i = 0; i < A.size(); ++i) {
        v[i] = A[i][n];
    }
    return v;
}
Intor Column_of_Intrix(Intrix A, int n)
{
    Intor v(A.size());
    for (int i = 0; i < A.size(); ++i) {
        v[i] = A[i][n];
    }
    return v;
}
void Save_Vector(const string& fn, const Vector& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(double e : v)   fil << e << endl;
}
void Save_Intor(const string& fn, const Intor& v)
{
    ofstream fil(fn);
    if(!fil)    { cout << "Filåbning mislykkedes."; return;   }
    for(int e : v)   fil << e << endl;
}
/*
intrix S_SS = Load_Dates("permnos_dates.csv");
Load SP500 //TODO
vector<double> beta_period(3);
for (int i = 0; i < 3; ++i)
{
    beta_period[i] = CovSP500(SP500, SP500, period_nr[0+i], period_nr[1+i]);
}
*/
void day_est_checker(int max)
{
    Intrix Dates = Load_Intrix("Dates_DR_Compressed.txt", max);
    Intor SP500_Dates = Load_Intor("List_of_dates.txt");

    if(max < 0 or max > Dates.size()) max = Dates.size();

    int True_date, Est_date;
    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][1];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][1], SP500_Dates)];
        if(True_date < Est_date)    cout << True_date << " <1 " << Est_date << ",  i = " << i << endl;
    }

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][2];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][2], SP500_Dates)];
        if(True_date < Est_date)    cout << True_date << " <2 " << Est_date << ",  i = " << i << endl;
    }

    cout << " --------------------------------------------- " << endl << endl;

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][1];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][1], SP500_Dates)];
        if(True_date > Est_date)    cout << True_date << " >1 " << Est_date << ",  i = " << i << endl;
    }

    for (int i = 0; i < max; ++i)
    {
        True_date = Dates[i][2];
        Est_date = SP500_Dates[Find_date_integer(Dates[i][2], SP500_Dates)];
        if(True_date > Est_date)    cout << True_date << " >2 " << Est_date << ",  i = " << i << endl;
    }
}

/*
Matrix Load_DR_Compressed(string fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Matrix(0);   }
    Matrix DR(36146);
    int n;

    for (int i = 0; i < 36146; ++i)
    {
        fil >> n;
        DR[i].resize(n);
        for (int j = 0; j < n; ++j) { fil >> DR[i][j]; }
    }
    return DR;
}
*/  // not using || pushback
/*
Matrix Load_DR_Compressed(string fn)
{
    ifstream fil(fn);
    if(!fil) {  cout << "Could not read the file " << fn << ".";  return Matrix(0);   }

    Matrix DR;
    double number;
    Vector Stock;
    int i=0;

    fil >> number;
    Stock.push_back(number);

    while(fil >> number)
    {
        if(number < 10000)  Stock.push_back(number);
        else
        {
            DR.push_back(Stock);
            Stock = Vector(0);
            Stock.push_back(number);
            if(i % 500 == 0)  cout << "ID = " << number << " i = " << i << endl;
            i++;
        }
    }
    return DR;
}
*/  // stock and vector || pushback