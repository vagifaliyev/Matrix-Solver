#pragma once

#include "Matrix.h"


template <class T>
class Symmetric: public Matrix<T>
{
// This class stores the top half of the symmetric matrix only
public:
    // Constructor
    Symmetric(int rows, int cols, bool preallocate);
    // Constructor where user preallocates memory
    Symmetric(int rows, int cols, T* values_ptr);
    // destructor
    virtual ~Symmetric();

    // turns dense matrix into symmetric format
    // "this" is the input, symm matrix to write to 
    // input variable is the dense matrix to be transformed 
    virtual void dense2symm(const Matrix<T> &sdense_mat);

    // turns symmetric matrix into dense format
    // "this" is the input, symm matrix to be transformed
    // input variable is the dense matrix to write to
    virtual void symm2dense(Matrix<T> &dense_mat);

    // Prints matrix in a symmetric fashion.
    // The user expierience is the same as for 
    // the parent class.
    virtual void printMatrix();


    // The memory allocated to b is assumed to be
    // of size to store as many values as the amount
    // columns of this matrix. This Matrix and b need
    // to be of same datatype.
    // The memory for x must be allocated BEFORE
    // calling this function.
    // A(this) * b = x
    virtual void matVecMult(T* &b, T* x);

    // matrix dimensions of this and mat_right
    // must be suitable for matrix multiplication 
    //void matMatMult(Symmetric<T> &mat_right, Symmetric<T> &output);
    
    // Linear Solver using SOR method
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
    // Where the residual res = b - A * x
    // or until reaching 10000 iterations
    // A(this)*x = b
    virtual void SOR(T *b, T *x, double omega, double atol);



  

    int len = -1;
};

