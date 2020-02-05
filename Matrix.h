#pragma once
#include <memory>

template<class T>
class Matrix
{
public:

    // constructor where we want to preallocate ourselves
    Matrix(int rows, int cols, bool preallocate);
    // constructor where we already have allocated memory outside
    Matrix(int rows, int cols, T* values_ptr);
    // destructor
    virtual ~Matrix();


    // FILL FUNCTIONS:

    // fills matrix values with random numbers
    // from 0 to 100
    virtual void fillRandom();
    // fills matrix values with zeros numbers
    virtual void fillZeros();
    // fills matrix values with ones numbers
    virtual void fillOnes();
    // Turns matrix into an identity matrix.
    // If a non-square matrix is input
    // the matrix will have values of one at
    // the diagonal starting at top-left corner of matrix.
    virtual void makeIdentity();

    // generates a Symmetric Positive Definite matrix
    // type = 1: dense matrix with high density filled with integers
    // type = 2: dense matrix with high density filled with doubles
    // type = 3: dense matrix where only the diagonal is populated
    // type = 4: dense matrix where only the tri-diagonal is populated
    // type = 5: dense matrix where only the penta-diagonal is populated
    // type = 6: dense matrix: diagonals populated: -5, -3, 0, 3, 5 
    //          (0: middle diag; negative: below middle diag; positive: above middle diag)
    virtual void generateSPD(int type);



    // PRINT FUNCTIONS

    // Print out the values in our matrix
    virtual void printValues();
    // Prints the matrix in a matrix-like fashion
    virtual void printMatrix();


    // PERFORM SOME OPERATIONS ON A MATRIX 

    // transpose function
    // the output matrix object needs to have rows and
    // columns switched with the this matrix
    virtual void transpose (Matrix<T>& x);
    // sums up all values inside matrix
    virtual void matSum(T &sum);
    // gets the trace of this matrix
    virtual void trace(T &trc);
    // counts the amount of non zero values
    virtual int countNonZeros();
    //computes the inverse of this matrix
    virtual bool inverse (Matrix<T> &inverse);
    // coFactor
    // used for changing the sign of alternate cells 
    // l and k represent the dimensions of the matrix
    // n is current dimension of M matrix
    virtual void coFactor (Matrix<T> &temp, int l, int k, int n); 
    //determinant
    // "this" represents the matrix that det is requested for
    // n is the current dimension of the matrox
    virtual double determinant (int n);
    //adjugate
    // "this" is the input matrix
    // variable is where to output the data to
    // swapping positions over the diagonal
    virtual void adjugate (Matrix<T> &adj);
    // compute the upper triangle reduced matrix given vector b
    // b is input vector
    // this is input matrix
    // out_mat is output matrix
    // x is output vector
    virtual void upper_triangle(T* b, Matrix<T> &out_mat, T* &x);
    // Back Substitution on the system Ax = b
    // Returns x, the solution
    // b is input vector
    // this is input matrix
    // x is output vector
    virtual void back_substitution(T* b, T* &x);
    // Forward Substitution on the system Ax = b
    // Returns x, the solution
    // b in input vector
    // this is input matrix
    // x is output vector
    virtual void for_substitution(T* b, T* &x);
    //LU decompisition 
    // "this" represents the input matrix
    //function will output lower and upper triangle matrices
    //which will be used by LUSolve
    virtual void luDecomposition(Matrix<T> &lower, Matrix<T> &upper);


    // VECTOR AND MATRIX MULTIPLICATIONS

    // The memory allocated to b is assumed to be
    // of size to store as many values as the amount
    // columns of this matrix. This Matrix and b need
    // to be of same datatype.
    // The memory for x must be allocated BEFORE
    // calling this funciton.
    // A(this) * b = x (x is the array to be populated by this function)
    virtual void matVecMult(T* &b, T* x);

    // matrix dimensions of this and mat_right
    // must be suitable for matrix multiplication
    // A(this) * B = X
    void matMatMult(Matrix<T> &B, Matrix<T> &X);
    // this function is not virtual because the imput for children classes
    // will be of different class

    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    ////////// LINEAR SOLVERS ////////////////////////////////

    // memory allocated to input arrays of all linear solvers
    // is asusmed to be preallocated by the user and to be the 
    // right length
    
    // linear solver by calculating the inverse of a matrix and 
    // then multiplying it with the vector b to get the solution x
    // A(this) * x = b
    virtual void linear_solver_inv(T* &b, T* x);

    // Linear Solver using Gaussian Elimination
    // A(this) * x = b
    virtual void gaussian_elim(T* b, T* &x);
    
    // Linear Solver using LU Decompisition
    // A(this) * x = b
    virtual void luSolve(T*b, T* &x);

   // Linear Solver Conjugate Gradient method
   // The function will iterate until reaching
   // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
   // Where the residual res = b - A * x
   // or until reaching 10000 iterations
   // A(this) * x = b
    virtual void conjugateGradient(T *b, T *x, double atol);
 
    // Linear Solver using SOR method
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
    // Where the residual res = b - A * x
    // or until reaching 10000 iterations
    // A(this)*x = b
    virtual void SOR(T *b, T *x, double omega, double atol);

    // Linear Solver using Chebyshev
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res1^2 + res2^2 + ...)
    // Where the residual res = b - A * x
    // or until reaching 10000 iterations
    // A(this)*x = b
    virtual void chebyshevIter(T* b, T* &x, T lMin, T lMax, double atol);


    // MULTIPLE LINEAR SOLVER
    // solve for multiple vectors b that represent the
    // vector space B. B is a matrix of size this->rows*m
    // where m is any natural number. The solution vector space X
    // is has the same dimensions as B.

    // Multi Linear Solver using Conjugate Gradient method
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res11^2 + res12^2 + ...)
    // Where the residual res = B - A * X
    // or until reaching 10000 iterations
    // A(this) * X = B
    void multiConjugateGradient(const Matrix<T> *B, Matrix<T> *X,\
    const double atol);
    // this function is not virtual because the imput for children classes
    // will be of different class

    // Multi Linear Solver using SOR method
    // The function will iterate until reaching
    // the absolut tolerance specified = sqrt(res11^2 + res12^2 + ...)
    // Where the residual res = B - A * X
    // or until reaching 10000 iterations
    // A(this) * X = B
    void multiSOR(const Matrix<T> *B, Matrix<T> *X, const double omega,\
    const double atol);
    // this function is not virtual because the imput for children classes
    // will be of different class


    // Explicitly using the C++11 nullptr here
    std::unique_ptr<T[]> values;
    int rows = -1;
    int cols = -1;

protected:
    bool preallocate = false;
// Private variables - there is no need for other classes 
// to know about these variables
private:
    int size_of_values = -1;
    
   
};
