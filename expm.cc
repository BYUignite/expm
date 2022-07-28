#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <vector>
#include <cmath>                   // log2, ceil, pow
#include <utility>                 // swap

using std::cout;
using std::endl;
using std::vector;

using boost::numeric::ublas::column_major;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::identity_matrix;

typedef matrix<double, column_major, std::vector<double> > bmatrix;

////////////////////////////////////////////////////////////////////////////////

extern "C" void dgesv_(int *n, int *nrhs, double *A, int *lda, int *ipiv,
                       double *b, int *ldb, int *info);

//--------- wrapper for dgesv_

void lsolve(bmatrix &A, bmatrix &BX){
    int n = A.size1();
    int nrhs = n;
    vector<int> ipiv(n);
    int info;
    dgesv_(&n, &nrhs, &A(0,0), &n, &ipiv[0], &BX(0,0), &n, &info);
    if(info != 0)
       std::cerr << std::endl << "ERROR: lsolve had info != 0: info = " << info << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void dgebal_(char* JOB, int *n, double *A, 
                        int *lda, int *ilo, int *ihi, 
                        double *scale, int *info);

//--------- wrapper for dgebal_

void balanceMatrix(bmatrix &A, vector<double> &scale, int &ilo, int &ihi) {
    int n = A.size1();
    int info;
    char JOB = 'B';           // B=Both scale and permute
    dgebal_(&JOB, &n, &A(0,0), &n, &ilo, &ihi, &scale[0], &info );
    //----------- fortran to c indexing
    ilo--;
    ihi--;
    for(int i=0; i<ilo; i++)
        scale[i]--;
    for(int i=ihi+1; i<n; i++)
        scale[i]--;

    if(info != 0)
       std::cerr << std::endl << "ERROR: balanceMatrix had info != 0: info = " << info << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void unBalanceMatrix(bmatrix &E, vector<double> &scale, int &ilo, int &ihi) {

    int n = E.size1();

    //------------------ undo the scaling

    double scj;
    for(int j=ilo; j<=ihi; j++) {
        scj = scale[j];
        for(int i=0; i<n; i++) {
            if(i==j) continue;
            E(j,i) *= scj;
            E(i,j) /= scj;
        }
    }

    //------------------ undo the permutations: reverse order of dgebal_

    int jj;
    if(ilo > 0)
        for(int j=ilo-1; j>=0; j--) {          // count down from ilo-1 to 0
            jj = int(scale[j]);                // scale[j] is the permutation index
            for(int k=0; k<n; k++)
                std::swap(E(k,j), E(k,jj));
            for(int k=0; k<n; k++)
                std::swap(E(j,k), E(jj,k));
        }
    if(ihi < n-1)
        for(int j=ihi+1; j<n; j++) {            // count up from ihi+1 to n-1
            jj = int(scale[j]);                 // scale[j] is the permutation index
            for(int k=0; k<n; k++)
                std::swap(E(k,j), E(k,jj));
            for(int k=0; k<n; k++)
                std::swap(E(j,k), E(jj,k));
        }
}

////////////////////////////////////////////////////////////////////////////////
// Exponential matrix
// Based on Julia's expm
// https://github.com/JuliaLang/julia/blob/d386e40c17d43b79fc89d3e579fc04547241787c/base/linalg/dense.jl#L395-L422
// Uses N. Higham "The scaling and squaring method of the matrix exponential revisited," SIAM Review Vol. 51, No. 4, pp. 747-764, 2009
// http://dx.doi.org/10.1137/090768539
// Tested against Matlab and Python expm
// Requires Boost and LAPACK
// AE: input/output: is A on input, is exp(A) on output
// D.O. Lignell July 27, 2022
// Compile Mac: g++ -std=c++11 -I/opt/homebrew/include -framework Accelerate -c expm.cc 

void expm(bmatrix &AE) {

    bmatrix &A = AE;
    bmatrix &E = AE;

    int n = A.size1();

    int ilo, ihi;
    vector<double> scale(n);
    balanceMatrix(A, scale, ilo, ihi);

    double nrmA   = norm_1(A);

    //-------------------------------------------------------------------------

    vector<double> C;

    //------------ For small nrmA, use lower order Padé-Approximations

    if(nrmA <= 2.1) {
        if(nrmA > 0.95)
            C = vector<double>{17643225600., 8821612800., 2075673600., 302702400.,
                                  30270240.,    2162160.,     110880.,      3960.,
                                        90.,          1.};
        else if(nrmA > 0.25)
            C = vector<double>{17297280.,8648640.,1995840.,277200.,
                                  25200.,   1512.,     56.,     1.};
        else if(nrmA > 0.015)
            C = vector<double>{30240.,15120.,3360.,
                                 420.,   30.,   1.};
        else
            C = vector<double>{120.,60.,12.,1.};

        bmatrix A2 = prod(A,A);

        bmatrix P(n,n,0.0);
        bmatrix U(n,n,0.0);
        bmatrix V(n,n,0.0);
        for(int i=0; i<n; i++) {
            P(i,i) = 1.0;
            U(i,i) = C[1];
            V(i,i) = C[0];
        }

        int K = int(C.size()/2) - 1;
        for(int k=1, k2=2; k<=K; k++, k2=2*k) {
            P = prod(P,A2);
            U += C[k2+1] * P;
            V += C[k2]   * P;
        }
        U = prod(A,U);
        P = V-U;
        E = U+V;

        lsolve(P, E);         // solving PE[:,i] = E[:,i] for each i
    }

    //------------ For large nrmA

    else {

        int si = std::ceil(std::log2(nrmA/5.4));   // dividing by power of 2 (2^s); undo below
        if(si > 0)
            A /= std::pow(2.0, si);

        C = vector<double>{64764752532480000., 32382376266240000., 7771770303897600.,
                            1187353796428800.,   129060195264000.,   10559470521600.,
                                670442572800.,       33522128640.,       1323241920.,
                                    40840800.,            960960.,            16380.,
                                         182.,                 1.};

        bmatrix A2 = prod(A,A);
        bmatrix A4 = prod(A2,A2);
        bmatrix A6 = prod(A2,A4);
        identity_matrix<double> I(n);

        bmatrix U = prod(A, prod(A6, C[13]*A6 + C[11]*A4 + C[9]*A2) + 
                            C[7]*A6 + C[5]*A4 + C[3]*A2 + C[1]*I);
        bmatrix V = prod(A6, C[12]*A6 + C[10]*A4 + C[8]*A2) + 
                    C[6]*A6 + C[4]*A4 + C[2]*A2 + C[0]*I;

        bmatrix P = V-U;
        E = V+U;

        lsolve(P, E);         // solving PE[:,i] = E[:,i] for each i

        if(si > 0)             // squaring (undo above divide by power of 2)
            for(int t=1; t<=si; t++)
                E = prod(E,E);
    }

    unBalanceMatrix(E, scale, ilo, ihi);
}
