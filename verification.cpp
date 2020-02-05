#include <iostream>
#include <ctime>
#include <memory>
#include <cmath>
#include <vector>
#include <iomanip>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "BandMatrix.h"
#include "BandMatrix.cpp"
#include "Symmetric.h"
#include "Symmetric.cpp"

using namespace std;

// parameters
int rows = 32;
int cols = 32;
double tolarance = 10e-10;

//compare arrays 
vector<int> compareArrays(double* a, double* b) 
{
    int n = rows;
    vector<int> diff;
    for(int i = 0; i < n; i++) 
    {
        if(abs(a[i]-b[i]) >= tolarance) 
        {
         diff.push_back(i);
        }
    }
    return diff;
}

//print vector
void printV(vector<int> input)
{
    for (int i = 0; i < input.size(); i++){
        cout << "#" << input[i] << " ";
    }
    // if empty means no difference
    if (input.size() == 0){
        cout << " all is within tolarance";
    }
}


int main()
{

    // initiate test matrices
    auto* test_mat1 = new Matrix<double>(rows, cols, true);
    auto* test_mat2 = new Matrix<double>(rows, cols, true);
    auto* test_mat3 = new Matrix<double>(rows, cols, true);
    auto* test_mat4 = new Matrix<double>(rows, cols, true);
    auto* test_mat5 = new Matrix<double>(rows, cols, true);
    auto* test_mat6 = new Matrix<double>(rows, cols, true);

    // generate test-SPDs for testing
    test_mat1->generateSPD(1);
    test_mat2->generateSPD(2);
    test_mat3->generateSPD(3);
    test_mat4->generateSPD(4);
    test_mat5->generateSPD(5);
    test_mat6->generateSPD(6);
    
    // cout << "type 1:" << endl;
    // test_mat1->printMatrix();
    // cout << "type 2:" << endl;
    // test_mat2->printMatrix();
    // cout << "type 3:" << endl; // uncomment to see test matrices
    // test_mat3->printMatrix();
    // cout << "type 4:" << endl;
    // test_mat4->printMatrix();
    // cout << "type 5:" << endl;
    // test_mat5->printMatrix();
    // cout << "type 6:" << endl;
    // test_mat6->printMatrix();

    //-------Test all DENSE mat solver--------\\

    // vector b    
    double* b1 = new double[rows];
    for (int i = 0; i < rows; i++)
    {
        b1[i] = rand() % rows + 0;
    }

    // INVERSE SOLVER ---- IS GREAT
    double* inv_out = new double[rows];
    clock_t start = clock();
    // test_mat1->linear_solver_inv(b1, inv_out);  // uncomment to see inverse solver
                                                   //- takes a while
    clock_t finish = clock();;
    double inv_time = ((double) (finish - start)) / CLOCKS_PER_SEC;
    
 
    // GAUSS ELIMINATION  ---- IS GREAT
    double* gauss_out = new double[rows];
    start = clock();
    test_mat1->gaussian_elim(b1, gauss_out);
    finish = clock();
    double gauss_time = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // LU SOLVE  --- IS GREAT
    double* lu_out = new double[rows];
    start = clock();
    test_mat1->luSolve(b1, lu_out);
    finish = clock();
    double lu_time = ((double) (finish - start)) / CLOCKS_PER_SEC;
    

    // CONJUGATE GRADIENT --- IS GREAT
    double* cg_out = new double[rows];
    start = clock();
    test_mat1->conjugateGradient(b1, cg_out, 10e-10);
    finish = clock();
    double cg_time = ((double) (finish - start)) / CLOCKS_PER_SEC;
    

    // SOR -----  IS GREAT
    double* sor_out = new double[rows];
    start = clock();
    test_mat1->SOR(b1, sor_out, 1, 10e-10);
    finish = clock();
    double sor_time = ((double) (finish - start)) / CLOCKS_PER_SEC;
    
    // CHEBYSHEV -----  IS GREAT
    double* cheb_out = new double[rows];
    start = clock();
    double emin1 = 1011.615;
    double emax1 = 2087.523;
    finish = clock();
    test_mat1->chebyshevIter(b1, cheb_out, emin1, emax1, 10e-10);
    double cheb_time = ((double) (finish - start)) / CLOCKS_PER_SEC;
    


    //-------Test all CRS mat solver--------\\


    // -----  IS GREAT
    int nnzs = test_mat6->countNonZeros();
    auto* test_sparse = new CSRMatrix<double>(rows, cols, nnzs, true);
    test_sparse->dense2csr(*test_mat6);

    // get dense result to compare
    double* cg_out_dense = new double[rows];
    start = clock();
    test_mat6->conjugateGradient(b1, cg_out_dense, 10e-10);
    finish = clock();
    double cg_time_dense = ((double) (finish - start)) / CLOCKS_PER_SEC;


    // CONJUGATE GRADIENT --- IS GREAT 
    double* cg_out_sparse = new double[rows];
    start = clock();
    test_sparse->conjugateGradient(b1, cg_out_sparse, 10e-10);
    finish = clock();
    double cg_time_sparse = ((double) (finish - start)) / CLOCKS_PER_SEC;
    
    // SOR ---- IS GREAT
    double* sor_out_sparse = new double[rows];
    start = clock();
    test_sparse->SOR(b1, sor_out_sparse, 1, 10e-10);
    finish = clock();
    double sor_time_sparse = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // CHEBYSHEV -- IS GREAT 
    double* cheb_out_sparse = new double[rows];
    start = clock();
    double emin2 = 597.730;
    double emax2 = 1526.657;
    finish = clock();
    test_sparse->chebyshevIter(b1, cheb_out_sparse, emin2, emax2, 10e-10);
    double cheb_time_sparse = ((double) (finish - start)) / CLOCKS_PER_SEC;

    //-------Test all BAND mat solver--------\\

    int bands = 3;
    auto* test_band = new BandMatrix<double>(rows, cols, bands, true);
    test_band->dense2band(bands, *test_mat4);

    // get dense result to compare
    double* cg_out_dense2 = new double[rows];
    start = clock();
    test_mat4->conjugateGradient(b1, cg_out_dense2, 10e-10);
    finish = clock();
    double cg_time_dense2 = ((double) (finish - start)) / CLOCKS_PER_SEC;


    // CONJUGATE GRADIENT -----  IS GREAT
    double* cg_out_band = new double[rows];
    start = clock();
    test_band->conjugateGradient(b1, cg_out_band, 10e-10);
    finish = clock();
    double cg_time_band = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // SOR -----  IS GREAT
    double* sor_out_band = new double[rows];
    start = clock();
    test_band->SOR(b1, sor_out_band, 1, 10e-10);
    finish = clock();
    double sor_time_band = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // CHEBYSHEV -- IS GREAT 
    double* cheb_out_band = new double[rows];
    start = clock();
    double emin3 = 602.774;
    double emax3 = 1522.364;
    finish = clock();
    test_band->chebyshevIter(b1, cheb_out_band, emin3, emax3, 10e-10);
    double cheb_time_band = ((double) (finish - start)) / CLOCKS_PER_SEC;



    //-------Test all SYMM mat solver--------\\


    auto* test_symm = new Symmetric<double>(rows, cols, true);
    test_symm->dense2symm(*test_mat4);

    // get dense result to compare
    double* cg_out_dense3 = new double[rows];
    start = clock();
    test_mat4->conjugateGradient(b1, cg_out_dense3, 10e-10);
    finish = clock();
    double cg_time_dense3 = ((double) (finish - start)) / CLOCKS_PER_SEC;


    // CONJUGATE GRADIENT -----  IS GREAT
    double* cg_out_symm = new double[rows];
    start = clock();
    test_symm->conjugateGradient(b1, cg_out_symm, 10e-10);
    finish = clock();
    double cg_time_symm = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // SOR -----  IS GREAT
    double* sor_out_symm = new double[rows];
    start = clock();
    test_symm->SOR(b1, sor_out_symm, 1, 10e-10);
    finish = clock();
    double sor_time_symm = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // CHEBYSHEV -- IS GREAT 
    double* cheb_out_symm = new double[rows];
    start = clock();
    double emin4 =  602.774;
    double emax4 = 1522.364;
    test_symm->chebyshevIter(b1, cheb_out_symm, emin4, emax4, 10e-10);
    finish = clock();
    double cheb_time_symm = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // ---- TEST MULTI DENSE SOLVERS ----\\


    // matrix b
    auto* b2 = new Matrix<double>(cols, 2, true);
    double* b3 = new double[rows];
    double* b4 = new double[rows];
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            b2->values[i * 2 + j] = rand() % rows + 0;
            if ((i * 2 + j) % 2 == 0)
            {
                b3[i] = b2->values[i * 2 + j];
            }
            else
            {
                b4[i] = b2->values[i * 2 + j];
            }
            
        }
    }

    //Conjugate gradient -----  IS GREAT
    // get dense result to compare
    //#x1
    double* cg_out_dense4 = new double[rows];
    start = clock();
    test_mat3->conjugateGradient(b3, cg_out_dense4, 10e-10);
    finish = clock();
    double cg_time_dense4 = ((double) (finish - start)) / CLOCKS_PER_SEC;

    // get dense result to compare
    //#x2
    double* cg_out_dense5 = new double[rows];
    start = clock();
    test_mat3->conjugateGradient(b4, cg_out_dense5, 10e-10);
    finish = clock();
    double cg_time_dense5 = ((double) (finish - start)) / CLOCKS_PER_SEC;

    //SOR -----  IS GREAT
    //#x1
    double* sor_out_dense4 = new double[rows];
    start = clock();
    test_mat3->SOR(b3, sor_out_dense4, 1, 10e-10);
    finish = clock();
    double sor_time_dense4 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    
    //#x2
    double* sor_out_dense5 = new double[rows];
    start = clock();
    test_mat3->SOR(b4, sor_out_dense5, 1, 10e-10);
    finish = clock();
    double sor_time_dense5 = ((double) (finish - start)) / CLOCKS_PER_SEC;


    // MULTI CONJUGATE GRADIENT -----  IS GREAT
    
    auto* cg_out_mult = new Matrix<double>(cols, 2, true);
    start = clock();
    test_mat3->multiConjugateGradient(b2, cg_out_mult, 10e-10);
    finish = clock();
    double cg_time_mult = ((double) (finish - start)) / CLOCKS_PER_SEC;
    
    // MULTI SOR -----  IS GREAT
    auto* sor_out_mult = new Matrix<double>(cols, 2, true);
    start = clock();
    test_mat3->multiSOR(b2, sor_out_mult, 1, 10e-10);
    finish = clock();
    double sor_time_mult = ((double) (finish - start)) / CLOCKS_PER_SEC;


    // PRINT DENSE SOLVER RESULTS
    cout << endl;
    cout << "----DENSE SOLVER RESULTS----" << endl;
    cout << endl << "no #" << "\t";
    cout << "Inverse" << "\t\t" << "Gauss" <<"\t\t" << "LU" << "\t\t";
    cout << "CG" << "\t\t" << "SOR" <<"\t\t"<< "Chebyshev" << endl;
    cout << endl;
    for (int i = 0 ; i < rows; i++)
    {
        cout << i << setw(16);
        cout << inv_out[i] << setw(16) << gauss_out[i] << setw(16);
        cout << lu_out[i] << setw(16)<< cg_out[i] << setw(16) << sor_out[i];
        cout << setw(16) << cheb_out[i] << endl;
    }
    cout << endl;
    cout << "Checked with tolarance of: " << tolarance << endl  << endl;
    cout << "Inverse solver took: " << inv_time << " s" << endl;
    cout << "check against CG: ";
    printV(compareArrays(inv_out, cg_out));
    cout << endl;
    cout << "Gauss elim solver took: " << gauss_time << " s" << endl;
    cout << "check against CG: ";
    printV(compareArrays(gauss_out, cg_out));
    cout << endl;
    cout << "L/U solver took: " << lu_time << " s" << endl;
    cout << "check against CG: ";
    printV(compareArrays(lu_out, cg_out));
    cout << endl;
    cout << "CG solver took: " << cg_time << " s" << endl;
    cout << "SOR solver took: " << sor_time << " s" << endl;
    cout << "check against CG: ";
    printV(compareArrays(sor_out, cg_out));
    cout << endl;
    cout << "Chebyshev solver took: " << cheb_time << " s" << endl;
    cout << "check against CG: ";
    printV(compareArrays(cheb_out, cg_out));
    cout << endl;




    // PRINT SPARSE SOLVER RESULTS
    cout << endl;
    cout << "----SPARSE SOLVER RESULTS----" << endl;
    cout << endl;
    cout << endl << "no #" << "\t";
    cout << "CG dense" << "\t" << "CG sparse" << "\t" << "SOR sparse";
    cout << "\t" << "Chebyshev sparse" <<  endl;
    cout << endl;
    for (int i = 0; i < rows; i++)
    {
        cout << i << setw(16) ;
        cout << cg_out_dense[i] << setw(16) << cg_out_sparse[i] << setw(16) ;
        cout << sor_out_sparse[i] << setw(16)  << cheb_out_sparse[i] <<  endl;
    }
    cout << endl;
    cout << "Checked with tolarance of: " << tolarance << endl;
    cout << "CG dense solver took: " << cg_time_dense << " s" << endl << endl;
    cout << "CG sparse solver took: " << cg_time_sparse << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cg_out_sparse, cg_out_dense));
    cout << endl;
    cout << "SOR sparse solver took: " << sor_time_sparse << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(sor_out_sparse, cg_out_dense));
    cout << endl;
    cout << "Chebyshev sparse solver took: " << cheb_time_sparse << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cheb_out_sparse, cg_out_dense));
    cout << endl;



    // PRINT BAND SOLVER RESULTS
    cout << endl;
    cout << "----BAND SOLVER RESULTS----" << endl;
    cout << endl;
    cout << endl << "no #" << "\t";
    cout << "CG dense" << "\t" << "CG band" << "\t\t" << "SOR band"<< "\t";
    cout << "Chebyshev band" << endl;
    cout << endl;
    for (int i = 0; i < rows; i++)
    {
        cout << i << setw(16);
        cout << cg_out_dense2[i] << setw(16) << cg_out_band[i] << setw(16);
        cout << sor_out_band[i] << setw(16) << cheb_out_band[i] <<  endl;
    }
    cout << endl;
    cout << "Checked with tolarance of: " << tolarance << endl;
    cout << "CG dense solver took: " << cg_time_dense2 << " s" << endl << endl;
    cout << "CG band solver took: " << cg_time_band << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cg_out_band, cg_out_dense2));
    cout << endl;
    cout << "SOR band solver took: " << sor_time_band << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(sor_out_band, cg_out_dense2));
    cout << endl;
    cout << "Chebyshev sparse solver took: " << cheb_time_band << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cheb_out_band, cg_out_dense2));
    cout << endl;



   // PRINT SYMM SOLVER RESULTS
    cout << endl;
    cout << "----SYMM SOLVER RESULTS----" << endl;
    cout << endl;
    cout << endl << "no #" << "\t";
    cout << "CG dense" << "\t" << "CG symm" << "\t\t" << "SOR symm";
    cout << "\t" << "Chebyshev symm" << endl;
    cout << endl;
    for (int i = 0; i < rows; i++)
    {
        cout << i << setw(16);
        cout << cg_out_dense3[i] << setw(16) << cg_out_symm[i] << setw(16);
        cout << sor_out_symm[i] << setw(16) << cheb_out_symm[i] <<  endl;
    }
    cout << endl;
    cout << "Checked with tolarance of: " << tolarance << endl;
    cout << "CG dense solver took: " << cg_time_dense2 << " s" << endl << endl;
    cout << "CG symm solver took: " << cg_time_symm << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cg_out_symm, cg_out_dense3));
    cout << endl;
    cout << "SOR symm solver took: " << sor_time_symm << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(sor_out_symm, cg_out_dense3));
    cout << endl;
    cout << "Chebyshev sparse solver took: " << cheb_time_symm << " s" << endl;
    cout << "check against CG_dense:";
    printV(compareArrays(cheb_out_symm, cg_out_dense3));
    cout << endl;

   // PRINT MULTI SOLVER RESULTS
    cout << endl;
    cout << "----MULTI SOLVER RESULTS----" << endl;
    cout << endl;
    cout << endl << "no #" << "\t";
    cout << "CG mult" << "\t\t\t\t" << "SOR mult" << endl;
    cout << "\t" << "x1" << "\t\t" << "x2";
    cout << "\t\t" << "x1" << "\t\t" << "x2";
    cout << endl;

    vector<double> cg_out_mult1;
    vector<double> cg_out_mult2;
    vector<double> sor_out_mult1;
    vector<double> sor_out_mult2;


    for (int i = 0; i < rows; i++)
    {
        cg_out_mult1.push_back(cg_out_mult->values[i * 2]);
        cg_out_mult2.push_back(cg_out_mult->values[i * 2 + 1] );
        sor_out_mult1.push_back(sor_out_mult->values[i * 2]);
        sor_out_mult2.push_back(sor_out_mult->values[i * 2 + 1]);
        cout << i << setw(16);
        cout << cg_out_mult1[i] << setw(16);
        cout << cg_out_mult2[i] << setw(16);
        cout << sor_out_mult1[i] << setw(16);
        cout << sor_out_mult2[i] << endl;


    }
    //turn vector into arrays
    double* cg_mult1= &cg_out_mult1[0];
    double* cg_mult2= &cg_out_mult2[0];
    double* sor_mult1= &sor_out_mult1[0];
    double* sor_mult2= &sor_out_mult2[0];

    cout << endl << endl << endl;
    cout << endl << "no #" << "\t";
    cout << "CG Dense" << "\t\t\t" << "SOR Dense" << endl;
    cout << "\t" << "x1" << "\t\t" << "x2";
    cout << "\t\t" << "x1" << "\t\t" << "x2";
    cout << endl << endl;
    
    for (int i = 0; i < rows; i++)
    {
        cout << i << setw(16);
        cout << cg_out_dense4[i] << setw(16);
        cout << cg_out_dense5[i] << setw(16);
        cout << sor_out_dense4[i] << setw(16);
        cout << sor_out_dense5[i] << endl;

    }

    cout << endl;
    cout << "Checked with tolarance of: " << tolarance << endl << endl;
    cout << "CG mult solver took: " << cg_time_mult << " s" << endl;
    cout << "check against CG_dense(x1):";
    printV(compareArrays(cg_mult1, cg_out_dense4));
    cout << endl;
    cout << "check against CG_dense(x2):";
    printV(compareArrays(cg_mult2, cg_out_dense5));
    cout << endl << endl;
    cout << "SOR mult solver took: " << sor_time_mult << " s" << endl;
    cout << "check against CG_dense(x1):";
    printV(compareArrays(sor_mult1, sor_out_dense4));
    cout << endl;
    cout << "check against CG_dense(x2):";
    printV(compareArrays(sor_mult2, sor_out_dense5));
    cout << endl;


    delete[] b1;
    delete b2;
    delete[] b3;
    delete[] b4;
    delete[] cheb_out_band;
    delete[] cheb_out_symm;
    delete[] sor_out_dense4;
    delete[] sor_out_dense5;
    delete cg_out_mult;
    delete sor_out_mult;

    delete[] inv_out;
    delete[] gauss_out;
    delete[] lu_out;
    delete[] cg_out;
    delete[] sor_out;
    delete[] cheb_out;
    delete[] cg_out_dense;
    delete[] cg_out_dense2;
    delete[] cg_out_dense3;
    delete[] cg_out_dense4;
    delete[] cg_out_dense5;
    delete[] cg_out_sparse;
    delete[] cg_out_band;
    delete[] cg_out_symm;
    delete[] cheb_out_sparse;
    delete[] sor_out_sparse;
    delete[] sor_out_band;
    delete[] sor_out_symm;

    delete test_mat1;
    delete test_mat2;
    delete test_mat3;
    delete test_mat4;
    delete test_mat5;
    delete test_mat6;

    delete test_sparse;
    delete test_band;
    delete test_symm;
}
