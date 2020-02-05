#pragma once

#include "Matrix.h"


template <class T>
class BandMatrix: public Matrix<T>
{
public:
    // constructor where we want to preallocate ourselves
    BandMatrix(int rows, int cols, int bands, bool preallocate);
    // constructor where we already have allocated memory outside
    BandMatrix(int rows, int cols, int bands, T *values_ptr);
    // destructor
    virtual ~BandMatrix();


    //Conversions
    // turns dense matrix into band format
    // "this" is the input, band matrix to write to 
    // input variable is the dense matrix to be transformed 
    virtual void dense2band(int bands, Matrix<T>& denseMat);

    // turns band matrix into dense format
    // "this" is the input, band matrix to be transformed
    // input variable is the dense matrix to write to
    virtual void band2dense(int bands, Matrix<T>& bandMat);
    
    // prints out the banded matrix
	virtual void printMatrix();


    // The memory allocated to b is assumed to be
    // of size to store as many values as the amount
    // columns of this matrix. This Matrix and b need
    // to be of same datatype.
    // The memory for x must be allocated BEFORE
    // calling this funciton.
    // A(this) * b = x
    virtual void matVecMult(T* &b, T* x);

/*
 * DID NOT HAVE ENOUGH TIME TO FULLY DEVELOP
 * WILL BE AVAILABLE NEXT PATCH
    // matrix dimensions of A(this) and B
    // must be suitable for matrix multiplication
    // Note: Unlike the Matrix matMatMult function, the object is created and 
    // assigned memory inside the function. This is because we do not know the
    // Number of Non-zero values before doing the multiplication.
    // PLEASE PASS ON EMPTY BandMatrix<T>* for X 
    // A(this) * B = X
    void matMatMult(BandMatrix<T> &B, BandMatrix<T>* &X);
    // this function is not virtual because the imput for children classes
    // will be of different class
*/

    // Linear Solver using SOR method
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
    // Where the residual res = b - A * x
    // or until reaching 10000 iterations
    // A(this)*x = b
    virtual void SOR(T *b, T *x, double omega, double atol);



    int bands = -1;
};

