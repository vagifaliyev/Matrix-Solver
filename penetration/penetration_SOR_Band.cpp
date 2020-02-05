#include <iostream>
#include <ctime>
#include <fstream>

#include "Matrix.h"
#include "Matrix.cpp"
#include "BandMatrix.cpp"
#include "BandMatrix.h"

using namespace std;

pair<int, double> penetration(int n, int type)
{
    int bands = 5;
    auto *dense_mat = new Matrix<double>(n, n, true);
    cout << "check for seg faults" << endl;
    dense_mat->generateSPD(5);
    cout << "generateSPD works" << endl;
    auto* penetration_mat = new BandMatrix<double>(n, n, bands, true);
    penetration_mat->dense2band(bands ,*dense_mat);
    cout << "dense2band works" << endl;
    double *b = new double[n];
    double *x = new double[n];
    for (int i = 0; i < n; i++)
    {
        b[i] = rand() % n + 0;
    }
    clock_t start = clock();
    penetration_mat->SOR(b, x, 1, 10e-10);
    clock_t finish = clock();
    double t = ((double) (finish - start)) / CLOCKS_PER_SEC;

    pair<int, double> result;
    result.first = n;
    result.second = t;

    delete dense_mat;
    delete penetration_mat;
    delete[] b;
    delete[] x;

    return result;
}

int main()
{
    int len = 6;
    int n[6] = {10, 100, 1000, 10000, 100000, 1000000};
    int type = 6;

    fstream data;
    data.open("./penetrationCGBand.dat");
    
    for (int i = 0; i < len; i++)
    {
        pair<int, double> data_point;
        data_point = penetration(n[i], type);
        data << data_point.first << data_point.second;
        cout << data_point.first << " " << data_point.second << endl;
        cerr << data_point.first << " " << data_point.second << endl;
    }
}
