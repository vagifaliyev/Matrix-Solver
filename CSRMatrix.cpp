#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "CSRMatrix.h"


using namespace std;

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate):\
Matrix<T>(rows, cols, false), nnzs(nnzs)
{
    // If we don't pass false in the initialisation list base constructor,
    // it would allocate values to be of size
    // rows * cols in our base matrix class
    // So then we need to set it to the real value we had passed in
    this->preallocate = preallocate;

    // If we want to handle memory ourselves
    if (this->preallocate)
    {
        // Must remember to delete this in the destructor
        this->values.reset(new T[this->nnzs]);
        this->row_position.reset(new int[this->rows+1]);
        this->col_index.reset(new int[this->nnzs]);
    }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr,\
int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr),\
nnzs(nnzs)
{
  this->row_position.reset(row_position);
  this->col_index.reset(col_index);
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{}

// turn a dense matrix intpo a CSR matrix
template<class T>
void CSRMatrix<T>::dense2csr(const Matrix<T>& denseMat)
{
    if (this->rows != denseMat.rows || this->cols != denseMat.cols)
    {
        cerr << "Dimensions of CSR matrix and dense Matrix don't match!" << endl;
        return;
    }
	this->row_position[0] = 0; // IA has the first element as 0
    int non_zero_count = 0;
    for (int i = 0; i < this->rows; i++)
    { 
		for (int j = 0; j < this->cols; j++)
        { 
			if (denseMat.values[i * this->rows + j] != 0)
            { 
                this->col_index[non_zero_count] = j;
                this->values[non_zero_count++] =\
                        denseMat.values[i * this->rows + j];
            }
		}
		row_position[i+1] = non_zero_count; 
	}
}

// turn a CSR matrix into a dense matrix
template<class T>
void CSRMatrix<T>::csr2dense(Matrix<T>& denseMat)
{
    if (this->rows != denseMat.rows || this->cols != denseMat.cols)
    {
        cerr << "Dimensions of CSR matrix and dense Matrix don't match!" << endl;
        return;
    }

    for (int i = 0; i < this->nnzs; i++)
    {
        for (int val_index = this->row_position[i];\
        val_index < this->row_position[i+1]; val_index++)
        {
            denseMat.values[i * this->cols + this->col_index[val_index]] =\
                    this->values[val_index];
        }

       
    }
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
    std::cout << "Printing CSR Matrix" << std::endl;
    std::cout << "Values: ";
    for (int j = 0; j< this->nnzs; j++)
    {  
        std::cout << this->values[j] << " ";      
    }
    std::cout << std::endl;
    std::cout << "row_position: ";
    for (int j = 0; j< this->rows+1; j++)
    {  
        std::cout << this->row_position[j] << " ";      
    }
    std::cout << std::endl;   
    std::cout << "col_index: ";
    for (int j = 0; j< this->nnzs; j++)
    {  
        std::cout << this->col_index[j] << " ";      
    }
    std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T* &b, T *x)
{
    if (b == nullptr || x == nullptr)
    {
        std::cerr << "b or x haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        x[i] = 0.0;
    }

    int val_counter = 0;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // Loop over all the entries in this col
        for (int val_index = this->row_position[i]; \
        val_index < this->row_position[i+1]; val_index++)
        {
            // This is an example of indirect addressing
            // Can make it harder for the compiler to vectorise!
            x[i] += this->values[val_index] * b[this->col_index[val_index]];

        }
    }

}



// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& B, CSRMatrix<T>* &X)
{
    // Check our dimensions match
    if (this->cols != B.rows)
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }
    
    // set temporary vectors to store indices and computed values
    vector<matRes<T>> temp;
    
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // Loop over the non-zero values in the i-th row of LHS matrix
        for (int val_index = this->row_position[i];\
        val_index < this->row_position[i+1]; val_index++)
        {
            // Loop over the non-zero values of RHS matrix
            // with matching indices:
            // row-index of LHS matrix matches column-index of RHS matrix
            for (int row_index = B.row_position[this->col_index[val_index]];\
            row_index < B.row_position[this->col_index[val_index] + 1];\
            row_index++)
            {
                // record position of the idividual values
                // ie record row index and column index for each multiplied value
                matRes<T> res; // create struct to store result
                // populate res object
                res.inner = B.col_index[row_index]; 
                res.outer = i;
                res.value = this->values[val_index] * B.values[row_index];
                // push to the temp struct
                temp.push_back(res);
            }
        }
    }


    // sort the entries of temp, initially 
    // according to column and then according to row
    sort(temp.begin(), temp.end(), this->compareMatResInner);
    sort(temp.begin(), temp.end(), this->compareMatResOut);
    // create vectors to temporarily store the arrays of our CSR output Matrix
    vector<T> values;
    vector<int> cols;
    vector<int> rows;
    rows.push_back(0); // first entry of rows is always 0
    // in this loop we essentially add together the results of our
    // multiplications to populate the output matrix
    for (int i = 0; i < temp.size(); i++)
    {
        // create temporary variable
        T tem = temp[i].value;
        // if the next value belongs to the same position in the output matrix
        // then add it, else move on.
        while ((i+1) != temp.size() && temp[i].outer == temp[i+1].outer &&\
        temp[i].inner == temp[i+1].inner)
        {
            tem += temp[i+1].value;
            i++;
        }
        // add it to the values array
        values.push_back(tem);
        // add the row position to the cols array
        cols.push_back(temp[i].inner);
        // get the number of rows that hold no entries
        if ((i+1) != temp.size() && (temp[i+1].outer - temp[i].outer > 0)) {
            int dif = temp[i+1].outer - temp[i].outer;
            // push to the rows the value size to indicate empty row
            for (int j = 0; j < dif; j++) rows.push_back(values.size());

        }
    }
    // make sure that rows with zeros in the end of 
    // the output matrix ger represented.
    while (rows.size() <= this->rows+1) {
        rows.push_back(values.size());
    }

    // initialize the pointer we gave as input
    X = new CSRMatrix<T>(this->rows, this->cols, values.size(), true);
    // populate the output matrix
    for (int i =0; i < values.size(); i++)
    {
        X->values[i] = values[i];
        X->col_index[i] = cols[i];
    }
    for (int i = 0; i < rows.size(); i++) X->row_position[i] = rows[i];
    
}

// Helper function for matMatMult
// compares the inner struct elements of two matRes objects
template<class T>
bool CSRMatrix<T>::compareMatResInner(matRes<T> a, matRes<T> b)
{
    // return bool for comparison
    return (a.inner < b.inner);
}

// Helper function for matMatMult
// compares the outer struct elements of two matRes objects
template<class T>
bool CSRMatrix<T>::compareMatResOut(matRes<T> a, matRes<T> b)
{
    // return bool for comparison
    return (a.outer < b.outer);
}


// the algorithm for this function is found on:
// https://en.wikipedia.org/wiki/Successive_over-relaxation
template<class T>
void CSRMatrix<T>::SOR(T *b, T *x, double omega, double atol)
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
        // res = b - Ax
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
            T diag_val = 0;
            // loop through all nnzs of the ith row
            for (int val_index = this->row_position[i];\
            val_index < this->row_position[i+1]; val_index++)
            {
                // if valie not in diagonal
                // multiply it with the element in x
                // that has the same index as the column index of the
                // non zero value of the current iteration
                if (this->col_index[val_index] != i)
                {
                    sigma += this->values[val_index] *\
                                x[this->col_index[val_index]];
                }
                else if(this->col_index[val_index] == i)
                {
                    diag_val = this->values[val_index];
                }

            }
            // computing ith value of result for current interation
            x[i] = (1 - omega) * x[i] +\
            (omega / diag_val) * (b[i] - sigma);
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
            cout << "SOR CRS is not converging!!!" << endl;
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
        } */
    }
    // print some cool stuff at the end
    cout << "Total iteration; " << count << endl;
    cout << "Total Residual: " << tot_rk << endl;

    delete[] x_A;
    delete[] r_k;

}

