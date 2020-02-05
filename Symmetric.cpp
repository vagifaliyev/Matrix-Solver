#include <iostream>
#include <vector>
#include <cmath>
#include "Symmetric.h"

using namespace std;

// Constructor - using an initialisation list here
template<class T>
Symmetric<T>::Symmetric(int rows, int cols, bool preallocate):\
Matrix<T>(rows, cols, false)
{
    if (this->rows != this->cols)
    {
        cout << "Only square matrices can be initiated as a Symmetric Matrix object" << endl;
        return;
    }
    this->preallocate = preallocate;
    this->len = this->rows * this->cols / 2 + (this->cols/2);
    if (preallocate)
    {
        this->values.reset(new T[len]);
    }
}

// Constructor - now just setting the value of our double pointer
template<class T>
Symmetric<T>::Symmetric(int rows, int cols, T *values_ptr):\
Matrix<T>(rows, cols, values_ptr)
{}

// destructor
template<class T>
Symmetric<T>::~Symmetric()
{}

template<class T>
void Symmetric<T>::dense2symm(const Matrix<T> &sdense_mat)
{
    vector<int> count;
    count.push_back(0);
    int idx;
    int nrows = this->rows;
    int ncols = this->cols;

    //assuming the user has checked that the dense matrix is already symmetric 
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols-i; j++) // only interested in above the diagonal 
        {
            this->values[count[i]+j] = sdense_mat.values[i*nrows+j+i]; //add according values to symm matrix
        }        
        idx = count[i] + this->cols - i; 
        count.push_back(idx); //update count index 
    }
}

template<class T>
void Symmetric<T>::symm2dense(Matrix<T> &dense_mat){
    vector<int> count;
    count.push_back(0);
    int idx;
    int nrows = this->rows;
    int ncols = this->cols;

    for (int i = 0; i < this->rows; i++)
    {
        for (int k = 0; k < i; k++)  //adding the reflections
        {
            dense_mat.values[i*nrows+k] = this->values[count[k]+i-k];

        }
        for (int j = 0; j < this->cols-i; j++) //adding the symm to a matrix 
        {
            dense_mat.values[i*nrows+j+i] = this->values[count[i]+j];
        }        
        idx = count[i] + this->cols - i;
        count.push_back(idx); // updating the index 
    }
}

// prints matrix in a symmetric fashion
template<class T>
void Symmetric<T>::printMatrix()
{
    // initiate a count to keep
    // to keep track of which index the 
    // next column starts with
    int count = 0;

    // loop through all rows
    for (int k = 0; k < this->rows; k++)
    {
        // fill all space with tabs until reaching
        // diaginal
        for (int j = 0; j < k; j++)
        {
            cout << "\t";
        }

        // fill the all spaces with the relevant entry
        // until reaching end of row
        for (int j = 0; j < this->cols - k; j++)
        {   
            // This is an example of referncing using a count
            cout << this->values[count + j] << "\t";
        }
        cout << endl;
        // update the count
        count += this->cols - k;
    }
}

// matVecMult
template<class T>
void Symmetric<T>::matVecMult(T* &vec_right, T* output)
{
    //ensure output vector is empty
    if (output != nullptr)
    {
        for (int i = 0; i < this->rows; i++)
        {
            output[i] = (T)0;
        }
    }

    vector<int> count; //keep track of the beggining index of each row 
    count.push_back(0); //initialize
    int idx; // pushing back idx to count 

    // looping through rows
    for (int i = 0; i < this->rows; i++)
    {
        // computing values left of diagonal
        for (int k = 0; k < i; k++)
        {
            output[i] += this->values[count[k]+i-k] * vec_right[k];
        }
        // computing values right of the diagonal
        for (int j = 0; j < this->cols-i; j++)
        {
            output[i] += this->values[count[i]+j] * vec_right[i+j];
        }
        // next index to push to count
        idx = count[i] + this->cols - i;
        count.push_back(idx);
    }
}


// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Successive_over-relaxation
template<class T>
void Symmetric<T>::SOR(T *b, T *x, double omega, double atol)
{
    // check if omega is set well
    if (omega <= 0 || omega >= 2)
    {
        cout << "Error: 0 < omega < 1" << endl;
        return;
    }
    // set default initial guess
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

    // initiate count
    int count = 0;
    // and the algorithm iterations start
    while (tot_rk >= atol)
    {   
        
        vector<int> ind_count;
        ind_count.push_back(0);
        // loop through rows and cols
        for (int i = 0; i < this->rows; i++)
        {
            double sigma = 0; // reset sigma to zero to re-compute it 
            for (int j = 1; j < this->cols - i; j++)
            {
                // recomputing sigma 
                // only adding terms outside the diagonal
                sigma += this->values[ind_count[i] + j] * x[i + j];
            }
            for (int j = 0; j < i; j++)
            {

                sigma += this->values[ind_count[j]+i-j] * x[j];

            }
            // computing ith value of result for current interation
            x[i] = (1 - omega) * x[i] +\
            (omega / this->values[ind_count[i]]) * (b[i] - sigma);

            // update count
            int idx = ind_count[i] + this->cols - i;
            ind_count.push_back(idx);
        }
        
        // reset residual
        tot_rk = 0;
        // compute A*x
        this->matVecMult(x, x_A);

        // compute total residual using vector norm        
        for (int i = 0; i < this->rows; i++)
        {
            r_k[i] = b[i] - x_A[i];
            tot_rk += r_k[i] * r_k[i];
        }
        tot_rk = sqrt(tot_rk); // sqrt to complete vector norm
        // update count
        count++;
        // break at maximum
        if (count == 10000)
        {
            cout << "It's not converging!!!" << endl;
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
