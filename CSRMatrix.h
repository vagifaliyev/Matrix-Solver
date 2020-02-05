#pragma once

#include "Matrix.h"
#include <memory>


// Helper struct for matMatMult
// it is used to store the multiplication results between 
// elements of two different matrices
// inner corresponds to column
// outer corresponds to row
// value corresponds to the result of the multiplication
template <typename T>
struct matRes
{
    int inner;
    int outer;
    T value;
};


template <class T>
class CSRMatrix: public Matrix<T>
{
public:

   // constructor where we want to preallocate ourselves
   CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
   // constructor where we already have allocated memory outside
   CSRMatrix(int rows, int cols, int nnzs, T *values_ptr,\
               int *row_position, int *col_index);
   // destructor
   virtual ~CSRMatrix();

   // Converts a dense matrix to a CSR matrix.
   // The function will be called with a common (dense) matrix object
   // and needs a CSR Matrix object as input.
   // The CSR matrix object needs to have
   // the same number of rows, cols and non-zero values as the imput matrix.
   // Please use the Matrix member function "countNonZeros"
   // to avoid coding unnecessarily.
   virtual void dense2csr(const Matrix<T>& denseMat);
   
   // Converts a CSR matrix to a dense matrix.
   // The function will be called with as a CSR matrix object
   // and needs a common (dense) Matrix object as input.
   // The common (dense) matrix object needs to have
   // the same number of rows and cols as this CSR matrix.
   virtual void csr2dense(Matrix<T>& denseMat);

   // Print out the values in our matrix
	virtual void printMatrix();


   // Perform some operations with our matrix
   
   // The memory allocated to b is assumed to be
   // of size to store as many values as the amount
   // columns of this matrix. This Matrix and b need
   // to be of same datatype.
   // The memory for x must be allocated BEFORE
   // calling this funciton.
   // A(this) * b = x
   virtual void matVecMult(T* &b, T *x);
   
   // matrix dimensions of A(this) and B
   // must be suitable for matrix multiplication
   // Note: Unlike the Matrix matMatMult function, the object is created and 
   // assigned memory inside the function. This is because we do not know the
   // Number of Non-zero values before doing the multiplication.
   // PLEASE PASS ON EMPTY CSRMatrix<T>* for X 
   // A(this) * B = X
   void matMatMult(CSRMatrix<T>& B, CSRMatrix<T>* &X);
   // this function is not virtual because the imput for children classes
   // will be of different class

   // Helper function for matMatMult
   // compares the inner struct elements of two matRes objects
   // usage described in: https://www.geeksforgeeks.org/sort-c-stl/
   static bool compareMatResInner(matRes<T> a, matRes<T> b);
   
   // Helper function for matMatMult
   // compares the outer struct elements of two matRes objects
   // usage described in: https://www.geeksforgeeks.org/sort-c-stl/
   static bool compareMatResOut(matRes<T> a, matRes<T> b);
   
   // Linear Solver using SOR method
   // The function will iterate until reaching
   // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
   // Where the residual res = b - A * x
   // or until reaching 10000 iterations
   // A(this)*x = b
   virtual void SOR(T *b, T *x, double omega, double atol);


   

   // Explicitly using the C++11 nullptr here
   std::unique_ptr<int[]> row_position;
   std::unique_ptr<int[]> col_index;
   

   // How many non-zero entries we have in the matrix
   int nnzs = -1;

// Private variables - there is no need for other classes 
// to know about these variables
private:
   
};

