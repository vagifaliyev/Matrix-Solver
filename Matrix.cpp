#include <iostream>
#include <typeinfo>
#include <cmath>
#include "Matrix.h"

using namespace std;

// Constructor - using an initialisation list here
template<class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate):\
rows(rows), cols(cols), size_of_values(rows * cols), preallocate(preallocate)
{
   // If we want to handle memory ourselves
   if (this->preallocate)
   {
      // Must remember to delete this in the destructor
      this->values.reset(new T[size_of_values]);
   }
}

// Constructor - now just setting the value of our double pointer
template<class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr):\
rows(rows), cols(cols), size_of_values(rows * cols)
{
    this->values.reset(values_ptr);
}

// destructor
template<class T>
Matrix<T>::~Matrix()
{}


// fills matrix values with random numbers
// between 0 and 100.
template<class T>
void Matrix<T>::fillRandom()
{
   for (int i = 0; i < this->rows * this->cols; i++)
   {
      this->values[i] = (T)(rand()%101 + 0);
   }
}

// fills matrix with zeros
template<class T>
void Matrix<T>::fillZeros()
{
   T val = (T)0;
   for (int i = 0; i < this->rows * this->cols; i++)
   {
      this->values[i] = val;
   }
}

// fills matrix with ones
template<class T>
void Matrix<T>::fillOnes()
{
   T val = (T)1;
   for (int i = 0; i < this->rows * this->cols; i++)
   {
      this->values[i] = val;
   }
}

// make an identity matrix
template<class T>
void Matrix<T>::makeIdentity()
{
   if (rows != cols)
   {
      cout << "A non-square matrix is made identity matrix. \
Matrix will be filled with zeros and have values of one \
at the diagonal starting at top-left corner of matrix." << endl;
   }
   T zero = (T)0;
   T one = (T)1;
   for (int i = 0; i < this->rows; i++)
   {
      for (int j = 0; j < this->cols; j++)
      {
         if (i == j)
         {
            this->values[i * this->cols + j] = one;
         }
         else
         {
            this->values[i * this->cols + j] = zero;
         }
         
      }
   }
}

template<class T>
void Matrix<T>::generateSPD(int type)
{
    // check if entry matrix is squared
    if (this->rows != this->cols)
    {
        cout << "Matrix needs to have the dimensions n x n" << endl;
        return;
    }

    // start switch statement to define which type if matrix to do
    switch (type)
    {
        case 1: // dense matrix with high density filled with integers
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i < j)
                    {   
                        // set random values for strictly upper triangle
                        this->values[i * this->cols + j] = this->rows - \
                                (rand() % this->rows + 0);
                    }
                    else if (j < i)
                    {
                        // set lower triangle as the upper on to achive symmetry
                        this->values[i * this->cols + j] = \
                                this->values[j * this->cols + i];
                    }
                    else
                    {
                        // set diagonal a random value guaranteed to be greater 
                        // than the sum of the rest of elements in ith row
                        this->values[i * this->cols + j] = this->rows * \
                                this->cols * 2 - \
                                (rand() % (this->rows*this->cols) + 0);
                    }
                }
            }
            break;
        }
        case 2: // dense matrix with high density filled with doubles
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i < j)
                    {
                        // set random double values for strictly upper triangle
                        double val =  (rand() % (this->rows * this->cols) + 0);
                                this->values[i * this->cols + j] = val / this->rows;
                    }
                    else if (j < i)
                    {
                        // set lower triangle as the upper on to achive symmetry
                        this->values[i * this->cols + j] = \
                                this->values[j * this->cols + i];
                    }
                    else
                    {
                        // set diagonal a random value guaranteed to be greater 
                        // than the sum of the rest of elements in ith row
                        this->values[i * this->cols + j] = \
                                this->rows * this->cols * 2 -\
                                (rand() % (this->rows*this->cols) + 0);
                    }
                }
            }
            break;
        }
        case 3: // dense matrix where only the diagonal is populated
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i == j)
                    {
                        // set diagonal values to random numbers
                        // (that though depend on the size of the matrix)
                        this->values[i * this->cols + j] =\
                                this->rows * this->cols * 1.5 -\
                                (rand() % (this->rows*this->cols) + 0);
                    }
                }
            }
            break;
        }
        case 4: // dense matrix where only the tri-diagonal is populated
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i == j)
                    {
                        // set diagonal values to random numbers
                        // (that though depend on the size of the matrix)
                        this->values[i * this->cols + j] =\
                                this->rows * this->cols * 1.5 -\
                                (rand() % (this->rows*this->cols) + 0);
                    }
                    else if (j - i == 1)
                    {
                        // set closest upper off diagonal to random values 
                        // between 0 and the amount of rows
                        // to guarantee SPD
                        this->values[i * this->cols + j] = this->rows -\
                                (rand() % (this->rows) + 0);
                    }
                    else if (i - j == 1)
                    {
                        // make symmetry by copying values over to the lower triangle
                        this->values[i * this->cols + j] =\
                                this->values[j * this->cols + i];
                    }
                }
            }
            break;
        }
        case 5:  // dense matrix where only the penta-diagonal is populated
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i == j)
                    {
                        // set diagonal values to random numbers
                        // (that though depend on the size of the matrix)
                        this->values[i * this->cols + j] =\
                                this->rows * this->cols * 1.5 -\
                                (rand() % (this->rows*this->cols) + 0);
                    }
                    else if (j - i == 1 || j - i == 2)
                    {
                        // set closest two upper off diagonals to random values 
                        // between 0 and the amount of rows
                        // to guarantee SPD
                        this->values[i * this->cols + j] = this->rows -\
                                (rand() % (this->rows) + 0);
                    }
                    else if (i - j == 1 || i - j == 2)
                    {                        
                        // make symmetry by copying values over to the lower triangle
                        this->values[i * this->cols + j] =\
                                this->values[j * this->cols + i];
                    }
                }
            }
            break;
        }
        case 6: // dense matrix: diagonals populated: -5, -3, 0, 3, 5
        {
            for (int i = 0; i < this->rows; i++)
            {
                for (int j = 0; j < this->cols; j++)
                {
                    if (i == j)
                    {
                        // set diagonal values to random numbers
                        // (that though depend on the size of the matrix)
                        this->values[i * this->cols + j] =\
                                this->rows * this->cols * 1.5 -\
                                (rand() % (this->rows*this->cols) + 0);
                    }
                    else if (j - i == 3 || j - i == 5)
                    {
                        // set closest third and fifth upper off diagonals to random values 
                        // between 0 and the amount of rows
                        // to guarantee SPD
                        this->values[i * this->cols + j] = this->rows -\
                                (rand() % (this->rows) + 0);
                    }
                    else if (i - j == 3 || i - j == 5)
                    {
                        // make symmetry by copying values over to the lower triangle
                        this->values[i * this->cols + j] =\
                                this->values[j * this->cols + i];
                    }
                }
            }
            break;
        }
        default:
            cout << "generateSPD: Type inserted does not exist!!!" << endl;
            break;
    }
}

// Just print out the values in our values array
template<class T>
void Matrix<T>::printValues() 
{ 
   std::cout << "Printing values" << std::endl;
	for (int i = 0; i< this->size_of_values; i++)
   {
      std::cout << this->values[i] << " ";
   }
   std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template<class T>
void Matrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   // looping through rows of matrix
   for (int j = 0; j< this->rows; j++)
   {  
      std::cout << std::endl;
      // looping through each element inside a row of the matrix
      for (int i = 0; i< this->cols; i++)
      {
         // We have explicitly used a row-major ordering here
         std::cout << this->values[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
}

// transpose function
template<class T>
void Matrix<T>::transpose (Matrix<T> &x)
{
    // check if intput matrix is squared 
    if (this->cols != x.rows || this->rows != x.cols)
    {
        cerr << "Dimensions don't match" << endl;
        return;
    }
    // switch the position of each element with its symmetric partner
    for( int i=0; i < this->rows; ++i)
    {
        for(int j=0; j < this->cols; ++j) 
        {
        x.values[j * x.cols + i] = this->values[i * x.rows + j];
        }
    }
}

//Makes the sum of the matrix
template<class T>
void Matrix<T>::matSum(T &sum)
{
    sum = (T)0;
    for (int i = 0; i < this->rows * this->cols; i++)
    {
        sum += this->values[i];
    }
}

// calculates the trace of a matrix
template<class T>
void Matrix<T>::trace(T &trc)
{  
    if (this->cols != this->rows)
    {
        cerr << "matrix must be square matrix" << endl;
        return;
    }

    trc = (T)0;
    cerr << "THIS IS HOW FAR I GOT" << endl;

    // create identity matrix
    unique_ptr<Matrix<T>> identity(new Matrix<T>(this->rows, this->cols, true));
    identity->makeIdentity();

    unique_ptr<Matrix<T>> trace_mat(new Matrix<T>(this->rows, this->cols, true));

    for (int i = 0; i < this->rows * this->cols; i++)
    {
        trace_mat->values[i] = identity->values[i] * this->values[i];
    }
    trace_mat->matSum(trc);
}

template<class T>
int Matrix<T>::countNonZeros()
{
    int nnzs = 0;
    for (int i = 0; i < this->rows * this->cols; i++)
    {
        if (this->values[i] != 0)
        {
            nnzs++;
        }
    }
    return nnzs;
}
  
/* Function to calculate and store inverse
with inverse(M) = adj(M)/det(M) formula
returns error if the matrix is singular */
template<class T>
bool Matrix<T>::inverse (Matrix<T> &inverse) 
{
    // check if matrix is squared 
    // if not throw error
    if (this->cols != this->rows)
    {
        cerr << "Only square matrices are invertible" << endl;
        return false;
    }
    int N = this->rows;
    // Find determinant of M[][] 
    double det = this->determinant(N);
    if (det == 0) 
    { 
        cout << "Inverse of a matrix exists only if the matrix is non-singular"; 
        return false; 
    } 
  
    // call Adjugate 
    unique_ptr<Matrix<T>> adj(new Matrix<T>(N, N, true)); 
    this->adjugate(*adj); 
  
    // Finding inverse using formula
    for (int i=0; i<N; i++) 
        for (int j=0; j<N; j++) 
            inverse.values[i* this->rows + j] = adj->values[i * adj->rows + j]/det;
    return true; 
} 

/* Applying a "checkerboard" of minuses to the "Matrix of Minors". 
In other words, we need to change the sign of alternate cells
Function to get cofactor of M in temp
where n is current dimension of M matrix */
template<class T>
void Matrix<T>::coFactor (Matrix<T> &temp, int l, int k, int n) 
{ 
    int i = 0, j = 0;
  
    // Looping thru each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != l && col != k) 
            { 
                temp.values[i * temp.rows + j] = \
                        this->values[row * this->rows + col];
                j++;
                // Row is filled so 
                if (j == n - 1) 
                { 
                    j = 0; //reset coloumn index
                    i++; // increase row index
                } 
            } 
        } 
    } 
} 
  
// Recursive function for finding determinant of matrix.
template<class T>
double Matrix<T>::determinant (int n)
{ 
    double D = 0; // Result initialised 
  
    // For 1x1 matrix, contains a single element 
    if (n == 1) 
        return this->values[0]; 
  
    unique_ptr<Matrix<T>> temp(new Matrix<T>(n, n, true)); // Store cofactors
  
    double sign = 1.0;  // Store sign multiplier 
  
    // Iterate for each element of first row 
    for (int z = 0; z < n; z++) 
    { 
        // Getting Cofactor of M[0][z] 
        this->coFactor(*temp, 0, z, n); 
        double det = temp->determinant(n-1);
        //temp->printMatrix();
        D += sign * this->values[z] * det; 
  
        // Alternating signs
        sign = -sign; 
    } 
    return D; 
} 
  
// Function to get Adjugate of M[N][N] in adj[N][N]. 
template<class T>
void Matrix<T>::adjugate (Matrix<T> &adj) 
{ 
    int N = this->rows;
    // For 1x1 matrix, contains a single element 
    if (N == 1) 
    { 
        adj.values[0] = 1; 
        return; 
    } 
  
    // temp is used to store cofactors of M[][] 
    int sign = 1;
    unique_ptr<Matrix<T>> temp(new Matrix<T>(N, N, true)); // Store cofactors 

    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
        { 
            // Get cofactor of M[i][j] 
            this->coFactor(*temp, i, j, N); 
  
            // If sum of row and column indexes is even,
            // sign of adj[j][i] is positive
            sign = ((i+j)%2==0)? 1: -1; 
 
            // transpose of the cofactor matrix, i.e interchanging rows and columns
            adj.values[j * this->cols + i] = (sign)*(temp->determinant(N-1)); 
        } 
    }
} 



// compute the upper triangle reduced matrix given vector b
template<class T>
void Matrix<T>::upper_triangle(T* b, Matrix<T> &out_mat, T* &x)
{
  // check if the matrix is a square matrix
  // must be square
  if (this->rows != this->cols) {
    std::cerr << "Matrix is not square!" << std::endl;
  }

  // set our output matrix to this matrix
  for (int i = 0; i < this->rows * this->cols; i++) 
  {
      out_mat.values[i] = this->values[i];
  }
  // set the output vector to the input matrix
  for (int i = 0; i < this->rows; i++) x[i] = b[i];
  
  double s = 0.0;  // define constant
  for (int k = 0; k < (this->rows)-1; k++)
  {
    for (int i = k+1; i < this->rows; i++)
    {
        // set constant value
        s = out_mat.values[i * this->rows + k] /\
                out_mat.values[k * this->rows + k];
        for (int j = k; j < this->rows; j++) 
        {
            // update output matrix
            out_mat.values[i * this->rows + j] =\
                    out_mat.values[i * this->rows + j] -\
                    s * out_mat.values[k * this->rows + j];
        }
        // update output vector
        x[i] = x[i] - s * x[k];
    }
  }
}

// Back Substitution on the system Ax = b
// Returns x, the solution
template<class T>
void Matrix<T>::back_substitution(T* b, T* &x)
{
  double s;
  // set output vector to 0
  for (int i = 0; i < this->rows; i++) x[i] = 0;

  // iterate over rows
  for (int i = (this->rows)-1; i > -1; i--)
  {
    // set constant
    s = 0.0;
    // iterate over columns
    for (int j = i+1; j < this->cols; j++)
    {
        // update constant value
        s += this->values[i * this->rows + j] * x[j];
    }
    // update the output matrix
    x[i] = (b[i] - s)/this->values[i * this->rows + i];
  }
}

// Forward Substitution on the system Ax = b
// Returns x, the solution
template<class T>
void Matrix<T>::for_substitution(T* b, T* &x)
{
  double s;
  // set output vector to 0
  for (int i = 0; i < this->rows; i++) x[i] = 0.0;

  // iterate over rows
  for (int i = 0; i < this->rows; i++) {
    s = 0.0; // set constant
    for (int j = 0; j < i; j++) {
      // update constant
      s += this->values[i * this->rows + j] * x[j];
    }
    // update output vector
    x[i] = (b[i] - s)/this->values[i * this->rows + i];
  }
}


// N×N matrix A, assume that an LU decomposition exists. 
//I am using Doolittle's method as it provides an alternative way to factor 
//A into an LU decomposition without going through the hassle of Gaussian Elimination.
//                 A=LU 
//where L is an n×n lower triangular matrix whose main diagonal consists of 1s 
//and where U is an n×n upper triangular matrix)
//Reference: http://mathonline.wikidot.com/the-algorithm-for-doolittle-s-method-for-lu-decompositions
template<class T>
void Matrix<T>::luDecomposition(Matrix<T> &lower, Matrix<T> &upper){
    
    int N = this->rows;

    //initiliase Matrix
    for (int i=0; i<N*N; i++){
        lower.values[i] = 0;
        upper.values[i] = 0;
    }

	for (int i = 0; i < N; i++) { 
		// Lower Triangular 
		for (int k = i; k < N; k++) { 
			if (i == k) 
				lower.values[i*N+i] = 1; // Diagonals as 1 
			else { 
				double sum = 0;
                // sum of L(k,j) * U(j, i) 
				for (int j = 0; j < i; j++) 
					sum += (lower.values[k*N+j] * upper.values[j*N+i]); 
				// evaluate lower
				lower.values[k*N+i] =\
                        (this->values[k*N+i] - sum) / upper.values[i*N+i]; 
			} 
         // Upper Triangular 
		for (int k = i; k < N; k++)
        { 
			double sum = 0;
            // sum of L(i,j) * U(j,k) 
			for (int j = 0; j < i; j++)
            {
				sum += (lower.values[i*N + j] * upper.values[j*N+k]); 
            }
			// evaluate upper
			upper.values[i*N + k] = this->values[i * N + k] - sum;
		} 
		} 
	} 
}
        

// Do matrix vector multiplication
// x = A(this) * b
template<class T>
void Matrix<T>::matVecMult(T* &b, T* x) 
{
   // initialize output vector to 0s
   for (int i = 0; i < this->rows; i++)
   {
      x[i] = 0;
   }
   // iterate over rows and coloumns
   for (int j = 0; j < this->rows; j++)
   {
      for (int i = 0; i < this->cols; i++)
      {
            x[j] += this->values[j * this->cols + i] * b[i];
      }
   }
}

// Do matrix matrix multiplication
// output = this * mat_right
template<class T>
void Matrix<T>::matMatMult(Matrix<T> &B, Matrix<T> &X)
{
   // Check our dimensions match
   if (this->cols != B.rows)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (X.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != X.rows || B.cols != X.cols)
      {
         cerr << "Input dimensions for output matrix doesn't match output" << endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      X.values.reset(new T[this->rows * B.cols]);
   }

   // Set values to zero before hand
   for (int i = 0; i < X.size_of_values; i++)
   {
      X.values[i] = 0;
   }

   for(int i = 0; i < this->rows; i++)
   {
      for(int k = 0; k < this->cols; k++)
      {
         for(int j = 0; j < B.cols; j++)
         {            
               X.values[i * X.cols + j] += this->values[i * this->cols + k] *\
                                            B.values[k * B.cols + j];
         }
      }
   }
}


/////////////////LINEAR SOLVERS\\\\\\\\\\\\\\\\\\


// linear solver by calculating the inverse of a matrix and
// then multiplying it with a vector to get the solution
template<class T>
void Matrix<T>::linear_solver_inv(T* &b, T* x)
{
  // create temporary matrix to store the inverse
  unique_ptr<Matrix<T>> temp(new Matrix<T>(this->rows, this->cols, true));
  // calculate the inverse
  this->inverse(*temp);
  // get the solution by multiplying the inverse with the input vector
  temp->matVecMult(b, x);
}

// Gaussian Elimination
template<class T>
void Matrix<T>::gaussian_elim(T* b, T* &x)
{

  // create a temporary matrix and vector to store the intermediate steps
  unique_ptr<Matrix<T>> temp_mat(new Matrix<T>(this->rows, this->cols, true));
  T* temp_vec = new T[this->rows];
  // get the upper triangle form of our input matrix
  this->upper_triangle(b, *temp_mat, temp_vec);
  // do back-substitution to get the solution vector
  temp_mat->back_substitution(temp_vec, x);
  delete[] temp_vec;
  return;
}

// LU decomposition solver
template<class T>
void Matrix<T>::luSolve(T* b, T* &x)
{
  unique_ptr<Matrix<T>> upper(new Matrix<T>(this->rows, this->cols, true));
  unique_ptr<Matrix<T>> lower(new Matrix<T>(this->rows, this->cols, true));
  T* temp_vec = new T[this->rows];

  this->luDecomposition(*lower, *upper);

  lower->for_substitution(b, temp_vec);

  upper->back_substitution(temp_vec, x);
  delete[] temp_vec;
  return;

}

// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
template<class T>
void Matrix<T>::conjugateGradient(T* b, T* x, double atol)
{
    // set default initial guess
    for (int i = 0; i < this->cols; i++)
    {
        x[i] = (T)1;
    }

    // declare some array pointers that we ll need in the algorithm
    T* r_k = new T[this->rows];
    T* x_A = new T[this->cols];
    T* p_k = new T[this->rows];

    // make the matVecMult to get the initial residual
    this->matVecMult(x, x_A);
    // declare initial residual
    double tot_rk = 0;
    // compute initial residual 
    for (int i = 0; i < this->rows; i++)
    {
        r_k[i] = b[i] - x_A[i]; // local residual
        tot_rk += r_k[i] * r_k[i]; // norm of residual
        p_k[i] = r_k[i]; // initiate p0
    }
    tot_rk = sqrt(tot_rk); // norm of initial residual

    // decclare alpha and beta
    // as variables in stack since only singular double
    double alpha_k;
    double beta_k;
    // initiate count
    int count = 0;
    // start iterating the algorithm
    while (tot_rk >= atol)
    {  
        // update counter
        count++;
        
        // initiate p_A as the product of pk and this->matrix
        double* p_A = new double[this->rows];
        this->matVecMult(p_k, p_A);

        // initiate array pointers denominator and counter 
        // for computing alpha
        double r_dot = 0;
        double p_Ap = 0;
        // get denominator and counter
        for (int i = 0; i < this->rows; i++)
        {
            r_dot += r_k[i] * r_k[i];
            p_Ap += p_k[i] * p_A[i];
        }
        // compute alpha
        alpha_k = r_dot / p_Ap;
                
        // getting x_k+1 and r_k+1
        // initiate array pointer to store rk+1
        double* r_k_nxt = new double[this->rows];
        // use the same loop to compute new result of kth iteration 
        // and rk+1
        for (int i = 0; i < this->rows; i++)
        {   
            x[i] += alpha_k * p_k[i];
            r_k_nxt[i] = r_k[i] - alpha_k * p_A[i];
        }

        // set residual variable to zero
        tot_rk = 0;
        // initialte denominator and counter 
        // to compute beta 
        double beta_top = 0;
        double beta_bot = 0;
        for (int i = 0; i < this->rows; i++)
        {
            // compute residual and beta rate in same loop
            tot_rk += pow(r_k_nxt[i], 2);
            beta_top += r_k_nxt[i] * r_k_nxt[i];
            beta_bot += r_k[i] * r_k[i];
        }
        // compu6te beta
        beta_k = beta_top / beta_bot;
        // take sqrt of residual to get the norm 
        tot_rk = sqrt(tot_rk);

        // compute pk and pass rk+1 to rk memory in same loop
        for (int i = 0; i < this->rows; i++)
        {
            p_k[i] = r_k_nxt[i] + beta_k * p_k[i];
            r_k[i] = r_k_nxt[i];
        }

        // print some nice stuff at the frequency you wish
        // change the integer after the % operator to change 
        // the frequency of printing
        /*
        if (count % 1 == 0)
        {
            cout << "Iteration: " << count << endl;
            cout << "Residual: " << tot_rk << endl;
            cout << "alpha: " << alpha_k << endl;
            cout << "beta: " << beta_k << endl;

            for (int i = 0; i < this->cols; i++)
            {
                cout << x[i] << endl;
            }
        }
        cout << endl;
        */

        delete[] p_A;
        delete[] r_k_nxt;
        
        // set limit of iterations for convergence
        if (count == 10000)
        {
            cout << "CG is not convergence,";
            cout << "be sure that the matrix is Symetric Positive Definite" << endl;
            break;
        }
    }
    // print some nice shit at the end
    cout << "Total iteration; " << count << endl;
    cout << "Total Residual: " << tot_rk << endl;

    delete[] r_k;
    delete[] x_A;
    delete[] p_k;
}



// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Successive_over-relaxation
template<class T>
void Matrix<T>::SOR(T *b, T *x, double omega, double atol)
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
        // loop through rows and cols
        for (int i = 0; i < this->rows; i++)
        {
            double sigma = 0; // reset sigma to zero to re-compute it 
            for (int j = 0; j < this->cols; j++)
            {
                // recomputing sigma 
                // only adding terms outside the diagonal
                if (i != j)
                {
                    sigma += this->values[i * this->cols + j] * x[j];
                }
            }
            // computing ith value of result for current interation
            x[i] = (1 - omega) * x[i] +\
            (omega / this->values[i * this->cols + i]) * (b[i] - sigma);
        }
        // reset residual
        tot_rk = 0;
        // compute A*x
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
            cout << "SOR is not converging!!!" << endl;
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


// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Chebyshev_iteration#cite_note-2
template<class T>
void Matrix<T>::chebyshevIter(T* b, T* &x, T lMin, T lMax, double atol)
{
    // initialize constants to be used in the algorithm
    double d = (lMax + lMin)/2;
    double c = (lMax - lMin)/2;

    // initialize our x and set it to our guess
    T* x0 = new T[this->cols];
    for (int i = 0; i < this->cols; i++)
    {
        x0[i] = (T)1;
    }

    // define output vector
    T* r = new T[this->cols];

    // r = b - A * x
    this->matVecMult(b, r);
    for (int i = 0; i < this->cols; i++)
    {
        r[i] = b[i] - r[i];
    }

    // initialize arrays for the iterative algorithm
    T* p = new T[this->cols];
    int i = 0;
    double alpha;
    double beta;
    double norm = atol * 2;
    while (norm >= atol)
    {
        if (i==0) {  // initial iteration
            // set p to r
            for (int j = 0; j < this->rows; j++) {
                p[j] = r[j];
            }
            alpha = 1.0/d; // define alpha
        } else if (i == 1) {  // second iteration
            beta = (1/2) * pow((c * alpha), 2); // beta = (1/2) * (c*alpha)^2
            alpha = 1/(d - beta/alpha);
            // p = r + beta * p
            for (int j = 0; j < this->rows; j++){
                p[j] = r[j] + beta * p[j];
            }
        } else { // the rest
            beta = pow(c*alpha/2, 2); // beta = (c*alpha/2)^2
            alpha = 1/(d - beta/alpha);
            // p = r + beta * p
            for (int j = 0; j < this->rows; j++) {
                p[j] = r[j] + beta*p[j];
            }
        }

        // x = x + alpha * p
        for (int j = 0; j < this->rows ; j++) {
            x0[j] = x0[j] + alpha * p[j];
        }

        // r = b - A*x;
        // calculate the norm of r
        this->matVecMult(x0, r);
        norm = 0.0;
        for (int j = 0; j < this->rows; j++) {
            r[j] = b[j] - r[j];
            norm += pow(r[j], 2);
        }
        norm = sqrt(norm);

        // max iterations check
        if (i == 10000) break;

        i++; // iterate
    }

    cout << "Total iteration; " << i << endl;
    cout << "Total Residual: " << norm << endl;

    // output result
    for (int j = 0; j < this->rows; j++) {
         x[j] = x0[j];
    }
    delete[] r;
    delete[] p;
    delete[] x0;

}


// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
template <class T>
void Matrix<T>::multiConjugateGradient(const Matrix<T> *B, Matrix<T> *X,\
const double atol)
{
    // check if dimensions of all matrices are right to 
    // be able to solve all linear systems
    if (this->cols != X->rows || this->rows != B->rows || X->cols != B->cols)
    {
        cout << "Dimensions are not right!!!" << endl;
        return;
    }

    // set default initial guess
    for (int i = 0; i < X->cols * X->rows; i++)
    {
        X->values[i] = (T)1;
    }

    // initiate some matrices we need during the algorithm 
    auto* r_k = new Matrix<T>(B->rows, B->cols, true);
    auto* p_k = new Matrix<T>(B->rows, B->cols, true);
    auto* x_A = new Matrix<T>(B->rows, B->cols, true);


    // precompute the matrix product A*x
    this->matMatMult(*X, *x_A);

    // initiate variable for residual
    // residual is the matrix norm 2
    double tot_rk = 0;
    for (int i = 0; i < B->rows * B->cols; i++)
    {
        r_k->values[i] = B->values[i] - x_A->values[i]; // compute local residual
        tot_rk += r_k->values[i] * r_k->values[i]; // get the square of local res
        p_k->values[i] = r_k->values[i]; // initiate p0
    }
    // get the sqrt to complete norm
    tot_rk = sqrt(tot_rk);

    // initiate arrays for alpah and beta
    double* alpha_k = new double[B->cols];
    double* beta_k = new double[B->cols];
    // initiate count 
    int count = 0;
    // start iterating the algorithm
    while (tot_rk >= atol)
    {  
        // update counter
        count++;
        
        // initiate matrix for A p matrix product and do the computation
        auto* p_A = new Matrix<T>(B->rows, B->cols, true);
        this->matMatMult(*p_k, *p_A);
        
        // initiate array pointers for dot product of rk and p*A*p
        double* r_dot = new double[r_k->cols];
        double* p_Ap = new double[p_k->cols];
        // set all allocated memory to zero
        for (int i = 0; i < B->cols; i++)
        {
            r_dot[i] = 0;
            p_Ap[i] = 0;
            alpha_k[i] = 0;
            beta_k[i] = 0;
        }

        // first getting denominator (dotproduct of rk)
        // and counter (p A p)
        for (int i = 0; i < B->rows; i++)
        {
            for (int j = 0; j < B->cols; j++)
            {
                r_dot[j] +=\
                        r_k->values[i * B->cols + j] * r_k->values[i * B->cols + j];
                p_Ap[j] +=\
                        p_k->values[i * B->cols + j] * p_A->values[i * B->cols + j];
            }
        }

        // then computing values for alpha
        // each value is assigned to one vector inside p_k matrix
        for (int j = 0; j < B->cols; j++)
        {
            alpha_k[j] = r_dot[j] / p_Ap[j];
        }

        // initiate matrix to store rk+1
        auto* r_k_nxt = new Matrix<T>(r_k->rows, r_k->cols, true);
        
        // compute the result for current iteration
        // and also the rk+1
        for (int i = 0; i < X->rows; i++)
        {   
            for (int j = 0; j < X->cols; j++)
            {
                X->values[i * X->cols + j] +=\
                        alpha_k[j] * p_k->values[i * X->cols + j];
                
                r_k_nxt->values[i * X->cols + j] =\
                        r_k->values[i * X->cols + j] -\
                        alpha_k[j] * p_A->values[i * X->cols + j];
            }
        }

        // set initial residual to zero to re-compute it
        tot_rk = 0;
        
        // initiate array pointers to store
        // denominator and counter to compute beta
        double* beta_top = new double[B->cols];
        double* beta_bot = new double[B->cols];
        // and set all allocated memory to zero
        for (int i = 0; i < B->cols; i++)
        {
            beta_top[i] = 0;
            beta_bot[i] = 0;
        }
        
        // making sure to take all the dot-product
        // of the individual vectors inside the
        // r_k / r_k+1 matrices
        for (int i = 0; i < B->rows; i++)
        {
            for (int j = 0; j < B->cols; j++)
            {
                tot_rk += pow(r_k_nxt->values[i * B->cols + j], 2);
             
                beta_top[j] += r_k_nxt->values[i * B->cols + j] *\
                        r_k_nxt->values[i * B->cols + j];

                beta_bot[j] += r_k->values[i * B->cols + j] *\
                        r_k->values[i * B->cols + j];
            }
        }
        tot_rk = sqrt(tot_rk); // sqrt to complete the matrix norm
        
        // compute beta with previously computed beta_top and beta_bot
        for (int j = 0; j < B->cols; j++)
        {
            beta_k[j] = beta_top[j] / beta_bot[j];
        }

        // compute p_k+1
        for (int i = 0; i < B->rows; i++)
        {
            for (int j = 0; j < B->cols; j++)
            {
                p_k->values[i * B->cols + j] =\
                        r_k_nxt->values[i * B->cols + j] +\
                        beta_k[j] * p_k->values[i * B->cols + j];
            }
        }

        // update r_k
        for (int i = 0; i < B->rows * B->cols; i++)
        {
            r_k->values[i] = r_k_nxt->values[i];
        }
        /*
        // print some nice stuff
        if (count % 1 == 0)
        {
            cout << endl << "Iteration: " << count << endl;
            cout << "Residual: " << tot_rk << endl;
            for (int i = 0; i < b.cols; i++)
            {   
                cout << "alpha, beta: ";
                cout << alpha_k[i] << " " << beta_k[i] << endl;
            }
            x.printMatrix();
            cout << endl;
        }*/
        
        delete[] r_dot;
        delete[] p_Ap;
        delete p_A;
        delete r_k_nxt;
        
        delete[] beta_bot;
        delete[] beta_top;
        
        // break if maximum number of iteration is reached
        if (count == 10000)
        {
            cout << "No convergence,";
            cout << "be sure that the matrix is Symetric Positive Definite" << endl;
            break;
        }
    }
    // print some more nice shit at the end
    cout << "Total iteration; " << count << endl;
    cout << "Total Residual: " << tot_rk << endl;

    delete x_A;
    delete r_k;
    delete p_k;
    delete[] alpha_k;
    delete[] beta_k;
}


// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Successive_over-relaxation
template<class T>
void Matrix<T>::multiSOR(const Matrix<T> *B, Matrix<T> *X, const double omega, const double atol)
{
    // check if omega is set well
    if (omega <= 0 || omega >= 2)
    {
        cout << "Error: 0 < omega < 1" << endl;
        return;
    }
    // check if dimensions of all matrices are right to 
    // be able to solve all linear systems
   if (this->cols != X->rows || this->rows != B->rows || X->cols != B->cols)
    {
        cout << "Dimensions are not right!!!" << endl;
        return;
    }

    // set default initial guess
    for (int i = 0; i < X->cols * X->rows; i++)
    {
        X->values[i] = (T)1;
    }

    // initiate Matrix pointers to store A*x
    // and rk
    auto* r_k = new Matrix<T>(B->rows, B->cols, true);
    auto* x_A = new Matrix<T>(B->rows, B->cols, true);

    // precompute the matrix product A*x
    this->matMatMult(*X, *x_A);

    // initiate variable for residual
    // residual is the matrix norm 2
    double tot_rk = 0;
    for (int i = 0; i < B->rows * B->cols; i++)
    {
        r_k->values[i] = B->values[i] - x_A->values[i]; // compute local residual
        tot_rk += r_k->values[i] * r_k->values[i]; // get the square of local res
    }
    // get the sqrt to complete norm
    tot_rk = sqrt(tot_rk);

    // initiate count
    int count = 0;
    // and the algorithm iterations start
    while (tot_rk >= atol)
    {
        // loop through rows and cols
        for (int i = 0; i < this->rows; i++)
        {
            double* sigma = new double[this->cols];
            // reset sigma to zero to re-compute it 
            for (int k = 0; k < X->cols; k++)
            {
                sigma[k] = 0;
            }
            for (int j = 0; j < this->cols; j++)
            {
                for (int k = 0; k < X->cols; k++)
                {
                    // recomputing sigma 
                    // only adding terms outside the diagonal
                    if (i != j)
                    {
                        
                        sigma[k] += this->values[i * this->cols + j] *\
                        X->values[i * this->cols + k];
                    }
                }
            }
            // computing ith value of result for current interation
            for (int k = 0; k < X->cols; k++)
            {
                X->values[i * X->cols + k] = \
                (1 - omega) * X->values[i * X->cols + k] +\
                (omega / this->values[i * this->cols + i]) *\
                (B->values[i * X->cols + k] - sigma[k]);
            }
            delete[] sigma;
        }
        // reset residual
        tot_rk = 0;
        // compute A*x
        this->matMatMult(*X, *x_A);

        // compute total residual using vector norm        
        for (int i = 0; i < B->rows * B->cols; i++)
        {
            // res = B - Ax
            // compute local residual
            r_k->values[i] = B->values[i] - x_A->values[i];
            // get the square of local res
            tot_rk += r_k->values[i] * r_k->values[i]; 
        }
        // get the sqrt to complete norm
        tot_rk = sqrt(tot_rk);

        // update count
        count++;
        // break at maximum
        if (count == 10000)
        {
            cout << "multiSOR is not converging!!!" << endl;
            break;
        }

     /*   // print some cool stuff
        cout << "Iteration: " << count << endl;
        cout << "Residual: " << tot_rk << endl;
        cout << "x:" << endl;
        for (int i = 0; i < x->rows; i++)
        {
            for (int j = 0; j < x->cols; j++)
            {
                cout << x->values[i * x->cols + j] << "  ";
            }
            cout << endl;
        }
        */

        
    }
    // print some cool shit at the end
    cout << "Total iteration; " << count << endl;
    cout << "Total Residual: " << tot_rk << endl;
    
    delete x_A;
    delete r_k;
}

