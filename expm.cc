#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <vector>
#include <cmath>                   // log2, ceil, pow, sqrt
#include <utility>                 // swap
#include <algorithm>               // max_element, max
#include <cfloat>                  // DBL_MAX

using std::cout;
using std::endl;
using std::vector;
using std::max_element;
using std::max;

using boost::numeric::ublas::column_major;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::identity_matrix;
using boost::numeric::ublas::row;

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

void balanceMatrix(bmatrix &A, vector<double> &scale, int &ilo, int &ihi, bool Lbalance) {

    int n = A.size1();

    if(Lbalance) {      // balance (scale)

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
    else {                // unbalance (unscale)

        //------------------ undo the scaling
        double scj;
        for(int j=ilo; j<=ihi; j++) {
            scj = scale[j];
            for(int i=0; i<n; i++) {
                if(i==j) continue;
                A(j,i) *= scj;
                A(i,j) /= scj;
            }
        }

        //------------------ undo the permutations: reverse order of dgebal_

        int jj;
        if(ilo > 0)
            for(int j=ilo-1; j>=0; j--) {          // count down from ilo-1 to 0
                jj = int(scale[j]);                // scale[j] is the permutation index
                for(int k=0; k<n; k++)
                    std::swap(A(k,j), A(k,jj));
                for(int k=0; k<n; k++)
                    std::swap(A(j,k), A(jj,k));
            }
        if(ihi < n-1)
            for(int j=ihi+1; j<n; j++) {            // count up from ihi+1 to n-1
                jj = int(scale[j]);                 // scale[j] is the permutation index
                for(int k=0; k<n; k++)
                    std::swap(A(k,j), A(k,jj));
                for(int k=0; k<n; k++)
                    std::swap(A(j,k), A(jj,k));
            }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Exponential matrix
// Based on Julia's expm (MIT license)
// https://github.com/JuliaLang/julia/blob/d386e40c17d43b79fc89d3e579fc04547241787c/base/linalg/dense.jl#L395-L422
// Nicholas J. Higham, "The squaring and scaling method for the matrix exponential revisited", 
// SIAM Journal on Matrix Analysis and Applications, 26(4), 2005, 1179-1193.
//      Reprinted online SIAM Review Vol. 51, No. 4, pp. 747-764, 2009 http://dx.doi.org/10.1137/090768539
//
// Tested against Matlab and Python expm
// Requires Boost and LAPACK
// AeA: input/output: is A on input, is exp(A) on output
// D.O. Lignell July 27, 2022
// Compile Mac: g++ -std=c++11 -I/opt/homebrew/include -framework Accelerate -c expm.cc 

void expm_higham05(bmatrix &AeA) {

    bmatrix &A  = AeA;
    bmatrix &eA = AeA;

    int n = A.size1();

    int ilo, ihi;
    vector<double> scale(n);
    balanceMatrix(A, scale, ilo, ihi, true);

    double nrmA   = norm_1(A);

    //-------------------------------------------------------------------------

    vector<double> b;

    //------------ For small nrmA, use lower order Padé-Approximations

    if(nrmA <= 2.1) {
        if(nrmA > 0.95)
            b = vector<double>{17643225600., 8821612800., 2075673600., 302702400.,
                                  30270240.,    2162160.,     110880.,      3960.,
                                        90.,          1.};
        else if(nrmA > 0.25)
            b = vector<double>{17297280.,8648640.,1995840.,277200.,
                                  25200.,   1512.,     56.,     1.};
        else if(nrmA > 0.015)
            b = vector<double>{30240.,15120.,3360.,
                                 420.,   30.,   1.};
        else
            b = vector<double>{120.,60.,12.,1.};

        bmatrix A2 = prod(A,A);

        bmatrix Q(n,n,0.0);
        bmatrix U(n,n,0.0);
        bmatrix V(n,n,0.0);
        for(int i=0; i<n; i++) {
            U(i,i) = b[1];
            V(i,i) = b[0];
        }

        int K = int(b.size()/2) - 1;
        for(int k=1, k2=2; k<=K; k++, k2=2*k) {
            Q = (k==1) ? A2 : prod(Q,A2);
            U += b[k2+1] * Q;
            V += b[k2]   * Q;
        }
        U  = prod(A,U);
        Q  = V-U;              // Higham's q
        eA = U+V;              // Higham's p here, then solved for r next, which we denote eA as in exp(A)

        lsolve(Q, eA);         // solving QeA[:,i] = eA[:,i] for each i
    }

    //------------ For large nrmA

    else {

        int si = std::ceil(std::log2(nrmA/5.4));   // dividing by power of 2 (2^s); undo below
        if(si > 0)
            A /= std::pow(2.0, si);

        b = vector<double>{64764752532480000., 32382376266240000., 7771770303897600.,
                            1187353796428800.,   129060195264000.,   10559470521600.,
                                670442572800.,       33522128640.,       1323241920.,
                                    40840800.,            960960.,            16380.,
                                         182.,                 1.};

        bmatrix A2 = prod(A,A);
        bmatrix A4 = prod(A2,A2);
        bmatrix A6 = prod(A2,A4);
        identity_matrix<double> I(n);

        bmatrix U = prod(A, prod(A6, b[13]*A6 + b[11]*A4 + b[9]*A2) + 
                            b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*I);
        bmatrix V = prod(A6, b[12]*A6 + b[10]*A4 + b[8]*A2) + 
                    b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*I;

        bmatrix Q = V-U;
        eA = V+U;

        lsolve(Q, eA);

        if(si > 0)             // squaring (undo above divide by power of 2)
            for(int t=1; t<=si; t++)
                eA = prod(eA,eA);
    }

    balanceMatrix(A, scale, ilo, ihi, false);
}

////////////////////////////////////////////////////////////////////////////////

// I. Najfeld and T.F. Havel, Derivatives of the Matrix Exponential and Their Computation, 
// Advances in Applied Mathematics, 16:321-375 (1995).
// This is referenced in Higham 2005
// AeA: input as A, output as exp(A)
// Summary:
// H(A) = Acoth(A) = A*(exp(2A)+I)/(exp(2A-I)) --> exp(2A) = (H(A)+A)/(H(A)-A)
// Let B = A/2^(d+1) = (A/2)/d^d --> 2B = A/2^d
// exp(2B) = (H(B)+B)/(H(B)-B)
// The diagonal Pade approximations of H are
// H(B) = M(B)/N(B)
//     H4(B) = (I + 4/9 * B^2 + 1/63 * B^4) / (I + 1/9 * B^2 + 1/945 * B^4)
// exp(2B) = (M/N + B)/(M/N - B) = (M + BN)/(M-BN)
// Solve (M-BN)exp(2B) = (M+BN) for exp(2B)
// Then square exp(2B) d times to get exp(A) = (exp(2B))^(2^d).
// Najfeld recommends H4 and a given gamma for single precision and H8 and a given gamma for double prec.
// In implementing this, minimize the number of intermediate matricies needed.
// doDP = true for double precision (default), false for single precision

void expm_najfeld_havel(bmatrix &AeA, bool doDP=true) {

    bmatrix &A    = AeA;               // using reference variables for notational convenience
    bmatrix &expA = AeA;

    const int n = A.size1();

    //----------- Balance matrix A

    int ilo, ihi;
    vector<double> scale(n);
    balanceMatrix(A, scale, ilo, ihi, true);

    //----------- Scale for Pade

    double  gamma4 = doDP ? 1.151922*4 : 0.9825211*4;
    bmatrix A2     = prod(A,A);
    double  nrmA2  = norm_1(A2);
    int     d      = (nrmA2 > gamma4) ? ceil(0.5*log2(nrmA2/gamma4)) : 0;

    //----------- Compute exp(2B)

    identity_matrix<double> I(n);

    bmatrix &B = AeA;
    B  = AeA * pow(2, -(d+1));
    bmatrix &B2 = A2;
    B2 = A2  * pow(2, -2*(d+1));
    bmatrix B4 = prod(B2,B2);

    bmatrix B6;
    bmatrix B8;
    if(doDP) {
        B6 = prod(B2, B4);
        B8 = prod(B4, B4);
    }

    bmatrix &BN = AeA;
    if(doDP) BN = prod(B, I + 7./51*B2 + 1./255*B4 + 2./69615*B6 + 1./34459425*B8);
    else     BN = prod(B, I + 1./9.*B2 + 1./945.*B4);

    bmatrix &M = B4;
    if(doDP) M  = I + 8./17*B2 + 7./255*B4 + 4./9945*B6 + 1./765765*B8;
    else     M  = I + 4./9.*B2 + 1./63. *B4;

    bmatrix &NN = A2;
    NN = M - BN;
    bmatrix &MM = AeA;
    MM = M + BN;

    lsolve(NN,MM);           // MM is now exp(2B), shares storage with expA, AeA

    //----------- Compute exp(A) by squaring d times

    for(int i=0; i<d; i++)                           // undo scaling
        expA = prod(expA, expA);

    //----------- un-balance matrix

    balanceMatrix(expA, scale, ilo, ihi, false);

}
