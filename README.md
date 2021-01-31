# ACSE 5: Matrix Library

This C++ library is meant for carrying out Matrix computations. Includes a total 4 types of matrices to use:

1. Dense Matrix:
saves strictly all values within an ordinary matrix.

2. CSR (sparse) Matrix:
consists of three arrays which save the (1) all non-zero values inside a matrix, (2) the row positioning and (3) the column indices of the non zero values.

3. Band Matrix:
saves only the values inside user-specified diagonals.

4. Symmetric Matrix:
saves only the upper triangle of of a matrix. 

All matrix types save their values in ROW MAJOR order.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

For the development of this library, only the standart c++11 libraries have been used.

```
#include <iostream>
#include <ctime>
#include <memory>
#include <cmath>
#include <vector>
#include <iomanip>
#include <vector>
#include <algorithm>
```

### Installation

- Clone this repo to your local machine using:

```
git clone https://github.com/vagifaliyev/Matrix-Solver
```

## Running the test

Test file runs all the avaliable linears solvers for the 4 matrix types.

Includes result comparison and time taken.

* Dense:
     * Inverse
     * Gauss
     * Lower–Upper (LU) decomposition 
     * Conjugate Gradient(CG) 
     * Successive over-relaxation (SOR)
     * Chebyshev
* Sparse(CSR):
     * Conjugate Gradient(CG) 
     * Successive over-relaxation (SOR)
     * Chebyshev
* Band:
     * Conjugate Gradient(CG) 
     * Successive over-relaxation (SOR)
     * Chebyshev
* Symmetric:
     * Conjugate Gradient(CG) 
     * Successive over-relaxation (SOR)
     * Chebyshev
* Multi-Solver(Dense only):
     * Conjugate Gradient(CG) 
     * Successive over-relaxation (SOR)

```
g++ Matrix.cpp CSRMatrix.cpp BandMatrix.cpp Symmetric.cpp verification.cpp -std=c++11
```

Note: The test call for inverse has been disabled for quicker testing as it is the most expensive function.

## Example Code

Header file contain in-depth explation on how to call the functions

### Creating a matrix 

```
auto* dense_mat1 = new Matrix<double>(cols, rows, true);

// fill in matrices and print them     
dense_mat1->fillRandom();
dens_mat1->printMatrix();


// this Matrix object has been implemented without
// preallocated memory. The int pointer matrix_data
// will pass on the memory in which the matrix values
// will be stored in.
auto* matrix_data = new double[rows * cols];
auto* dense_mat1 = new Matrix<double>(cols, rows, matrix_data);
  
```

### Call a function 

```
double* vec_right = new double[cols];

// set values of RHS vector in memory space
for (int i = 0; i < cols; i++)
{
    vec_right[i] = 2;
}

// initiate new poiter and allocate memory
// for result vector
double* vec_out = new double[rows];

// do matrix vector product
dense_mat->matVecMult(vec_right, vec_out);
```

## Documentation 

Please read [report.pdf](https://github.com/vagifaliyev/Matrix-Solver/blob/master/report.pdf) for details on our code of conduct, and algorthims.

Penetration analysis results can be viewed in folder [outputData](https://github.com/vagifaliyev/Matrix-Solver/tree/master/penetration/outputData)


## On the next Patch:

* Finish implementing MatMatMult for band & symmetric matrix 
* Parallize with openMP 
* Write a function to analyze the input matrix and direct it to the most aprropiate class.

## Acknowledgments

Thanks to the teaching staff for their enthusiastic classes and all the Teaching Assistants for their assistance.

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/vagifaliyev/Matrix-Solver/blob/master/LICENSE) file for details.
