#include <iostream>
#include <ctime>
#include <fstream>

#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

pair<int, double> penetration(int n, int type)
{
    auto *penetration_mat = new Matrix<double>(n, n, true);
    penetration_mat->generateSPD(2);
    double *b = new double[n];
    double *x = new double[n];
    for (int i = 0; i < n; i++)
    {
        b[i] = rand() % n + 0;
    }
    clock_t start = clock();
    penetration_mat->conjugateGradient(b, x, 10e-10);
    clock_t finish = clock();
    double t = ((double) (finish - start)) / CLOCKS_PER_SEC;

    pair<int, double> result;
    result.first = n;
    result.second = t;

    delete penetration_mat;
    delete[] b;
    delete[] x;

    return result;
}

int main()
{
    int len = 4;
    int n[4] = {10, 100, 1000, 10000};
    int type = 2;

    fstream data;
    data.open("./penetrationCG.dat");
    
    for (int i = 0; i < len; i++)
    {
        pair<int, double> data_point;
        data_point = penetration(n[i], type);
        data << data_point.first << data_point.second;
        cout << data_point.first << " " << data_point.second << endl;
    }
}