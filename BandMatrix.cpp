#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include "BandMatrix.h"

using namespace std;

// Constructor - using an initialisation list here
template <class T>
BandMatrix<T>::BandMatrix(int rows, int cols, int bands, bool preallocate):\
Matrix<T>(rows, cols, false), bands(bands)
{
    // If we don't pass false in the initialisation list base constructor,
    // it would allocate values to be of size
    // rows * cols in our base matrix class
    // So then we need to set it to the real value we had passed in
    this->preallocate = preallocate;

    if (this->cols != this->rows)
    {
        cerr << "matrix must be square matrix" << endl;
        return;
    }
    if (this->bands % 2 == 0)
    {
        cerr << "Only odd number of bands" << endl;
        return;
    }

    // If we want to handle memory ourselves
    if (this->preallocate)
    {
        // Must remember to delete this in the destructor
        this->values.reset(new T[this->rows * this->bands]);
        
    }
}

template <class T>
BandMatrix<T>::BandMatrix(int rows, int cols, int bands, T *values_ptr):\
Matrix<T>(rows, cols, values_ptr), bands(bands)
{}

// destructor
template <class T>
BandMatrix<T>::~BandMatrix()
{

    // The super destructor is called after we finish here
    // This will delete this->values if preallocate is true
}

template<class T>
void BandMatrix<T>::dense2band(int bands, Matrix<T>& denseMat) 
{
    
    int nrows = this->rows;
    int ncols = this->bands;
    // the bandwidth is usually given in odd number
    // hence adding 1 to get a real number
    int dx = (bands+1)/2; 

   for (int i = 0; i < nrows; i++) 
   {
        for (int j = 0; j < bands; j++) 
        {
            if (i==0) //initial border
            {
                if (j<dx-1)
                {
                    // add 0 for the missing elements 
                    this->values[i*nrows+j] = 0; 
                }
                else
                {
                    // add the rest of the elements
                    this->values[j] = denseMat.values[j+(1-dx)]; 
                }
            }
            else if (i==nrows-1) //final border 
            {
                if (j<=bands-dx)
                {
                    this->values[ncols*i+j] = denseMat.values[i*nrows+j+i+(1-dx)];
                }
            }
            else //middle domains
            {
                //move diagonally and add the elements according to bands
                this->values[ncols*i+j] = denseMat.values[i*nrows+j+i+(1-dx)]; 
            }
        }
    }
}

template<class T>
void BandMatrix<T>::band2dense(int bands, Matrix<T>& bandMat) 
{
    int nrows = this->rows;
    int ncols = this->cols;
    // the bandwidth is usually given in odd number
    // hence adding 1 to get a real number
    int dx = (bands+1)/2; 

    for (int i = 0; i < nrows; i++) 
    {
        for (int j = 0; j < ncols; j++) 
            {
                if (i==0) //initial border
                {
                    if (j<=bands-dx){
                        bandMat.values[j] = this->values[j+dx-1]; 
                    }
                    else 
                    {
                        // add 0 for out of bounds elements
                        bandMat.values[j] = 0.0;
                    }
                }
                else if (i==nrows-1) //final border 
                {
                    if (j>i-dx){
                        bandMat.values[i*nrows+j] = this->values[i*bands+j-i+dx-1];
                    }
                    else 
                    {
                        // add 0 for out of bounds elements
                        bandMat.values[i*nrows+j] = 0; 
                    }
                }
                else // middle/end domain 
                {
                    if (j>=i-1 && j<i+bands-1)
                    {
                        // add values diagonally according to bands
                        bandMat.values[i*nrows+j] = this->values[i*bands+j-i+dx-1]; 
                    }
                    else
                    {
                        // add 0 for out of bounds elements
                        bandMat.values[i*nrows+j+i-1] = 0; 
                    }
                }
            }
    }    
}

template <class T>
void BandMatrix<T>::printMatrix()
{
    cout << "Printing Banded Matrix:" << endl;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->bands; j++)
        {
            cout << this->values[i * this->bands + j] << " ";
        }
        cout << endl;
    }
}

template<class T>
void BandMatrix<T>::matVecMult(T* &b, T* x)
{
    // setting all output values to zero
    for (int i = 0; i < this->rows; i++) x[i] = 0;


    // initialize variables that hold positional
    // arguments
    int start;
    int end;
    int half = (this->bands-1)/2;
    for (int i = 0; i < this->rows; i++)
    {
        // set start and end
        start = max(0, i-half)-i+half+1;
        end = min(this->rows-1, i+half)-i+half+1;
        
        for (int j = start-1; j < end; j++)
        {
            // do multiplication
            x[i] += this->values[i * this->bands + j] * b[j+i-half];
        }
    }
}
/*
 * DID NOT HAVE ENOUGH TIME TO FULLY DEVELOP
 * WILL BE AVAILABLE NEXT PATCH
template <class T>
void BandMatrix<T>::matMatMult(BandMatrix<T> &B, BandMatrix<T>* &X)
{
    
    // checking if dimensions of all matrices are correct
    if (this->cols != B.rows)
    {
        cerr << "Dimensions of matrices don't match" << endl;
        return;
    }
    
    // getting the values for the midle column of the output
    // first and last excluded
    int half = (this->bands-1)/2;
    int start;
    int end;
    int N = b.rows;
    int M = b.bands;
    int i= M - 2;
    vector<vector<T>> out;
    cout << "half: " << half << endl;

    while(i >= 0) {
        vector<T> temp;
        int row = 0;
        int col = M - i - 1;
        for (int j=0; j < M-i; j++) {
             temp.push_back(b.values[row * M + col]);
            row++;
            col--;
        }
        for (int k = M-i; k < N; k++) temp.push_back(0);
        out.push_back(temp);
        i--;
    }

    i = N - 2;
    while (i >= 1) {
        int row = N - i - 1;
        int col = M - 1;
        vector<T> temp;
        for (int k = N-(N-i)+half; k < N; k++) temp.push_back(0);
        cout << min(i, M-1) << endl;
        int ind = min(i, M-1);
        for (int j = 0; j <= ind; j++) {
            temp.push_back(b.values[row * M + col]);
            row++;
            col--;
        }
        for (int k=0; k < i-M+1; k++) temp.push_back(0);
        out.push_back(temp);
        i--;
    }

    for (int i = 0; i < out.size(); i++) {
        for (int j = 0; j < out[i].size(); j++) {
            cout << out[i][j] << " ";
        }
        cout << endl;
    }

    T* final = new T[this->rows * this->cols];
    for (int i = 0; i < this->rows*this->cols; i++) {
        final[i] = 0;
    }
    for (int i = 0; i < this->rows; i++)
    {
        start = max(0, i-half)-i+half+1;
        end = min(this->rows-1, i+half)-i+half+1;

        for (int k = 0; k < out.size();k++)
        {
            for (int j = start-1; j < end; j++)
            {
                final[i*this->rows + k] += \
                        out[k][j+i-half] * this->values[i*this->bands + j];
            }
        }
    }

    int bands = 0;
    int a = this->cols/2;
    int d = 0;
    cout << "a: " << a << endl;
    cout << "cols: " << this->cols;
    while (d < this->cols) {
        if (final[a*this->cols + d] != 0.0) {bands++;}
        d++;
    }

    auto* denseForm = new Matrix<T>(this->rows, this->cols, final);
    denseForm->printMatrix();

   // cout << "hohohoh" << endl;
    //x = new BandMatrix<T>(this->rows, this->cols, bands, true);
   // cout << "hehehehehe" << endl;
    //x->dense2band(bands, *denseForm);
   // cout << "hahahahha" << endl;

}
*/

// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Successive_over-relaxation
template<class T>
void BandMatrix<T>::SOR(T *b, T *x, double omega, double atol)
{
    // check if omega is set well
    if (omega <= 0 || omega >= 2)
    {
        cout << "Error: 0 < omega < 1" << endl;
        return;
    }
    // set initial guess
    for (int i = 0; i < this->cols; i++)
    {
        x[i] = (T)1;
    }

    // initiate array pointers to store A*x
    // and rk
    T* x_A = new T[this->cols];
    T* r_k = new T[this->rows];
    // compute A*x
    this->matVecMult(x, x_A);
    // initiate a variable for the residual
    double tot_rk = 0;
    // compute residual as vector norm
    for (int i = 0; i < this->rows; i++)
    {
        r_k[i] = b[i] - x_A[i];
        tot_rk += r_k[i] * r_k[i];
    }
    tot_rk = sqrt(tot_rk); // sqrt to complete vector norm

    // define new array with extra slots at the beginning
    // and end for buffer
    unique_ptr<T[]> x_temp(new T[this->cols + this->bands - 1]);
    // set beginning buffer to zero
    for (int i = 0; i < (this->bands-1)/2; i++) x_temp[i] = 0;
    // set middle part to initial guess of x
    for (int i = (this->bands-1)/2; i < (this->bands-1)/2 + this->cols;\
    i++) x_temp[i] = x[i + (this->bands-1)/2];
    // set end buffer to zero
    for (int i = (this->bands-1)/2 + this->cols;\
    i < (this->bands-1) + this->cols; i++) x_temp[i] = 0;

    // initiate count
    int count = 0;
    // and the algorithm iterations start
    while (tot_rk >= atol)
    {
        // loop through rows and cols
        for (int i = 0; i < this->rows; i++)
        {
            double sigma = 0; // reset sigma to zero to re-compute it 
            for (int j = 0; j < this->bands; j++)
            {
                // recomputing sigma 
                // only adding terms outside the diagonal
                if (j != (this->bands-1)/2)
                {
                    // using the buffer for top and bottom 
                    // of the banded matrix
                    sigma += this->values[i * this->bands + j] * x_temp[i + j];
                }
            }
            // computing ith value of result for current iteration
            x_temp[i + (this->bands-1)/2] = (1 - omega) *\
                x_temp[i + (this->bands-1)/2] +\
                (omega / this->values[i * this->bands + (this->bands-1)/2]) *\
                (b[i] - sigma);
        }
        // reset residual
        tot_rk = 0;
        // copying the x_temp over to actual solution x
        for (int i = 0; i < this->cols; i++)
        {
            x[i] = x_temp[i + (this->bands-1)/2];
        }
        // compute A*x for residual
        this->matVecMult(x, x_A);

        // compute total residual using vector norm        
        for (int i = 0; i < this->rows; i++)
        {
            // res = b - Ax
            r_k[i] = b[i] - x_A[i];
            tot_rk += r_k[i] * r_k[i];
        }
        tot_rk = sqrt(tot_rk); // sqrt to complete vector norm
        // update count
        count++;
        // break at maximum
        if (count == 10000)
        {
            cout << "SOR Band is not converging!!!" << endl;
            break;
        }
        /*
        // print some cool stuff
        cout << "Iteration: " << count << endl;
        cout << "Residual: " << tot_rk << endl;
        cout << "x:" << endl;
        for (int i = 0; i < this->rows; i++)
        {
            cout << x[i] << endl;
        }*/
    }

    // print some cool stuff at the end
    cout << "Total iteration; " << count << endl;
    cout << "Total Residual: " << tot_rk << endl;

    delete[] x_A;
    delete[] r_k;
}
