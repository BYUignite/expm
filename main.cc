//g++ -std=c++11 main.cc

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>

using namespace boost::numeric::ublas;
using std::cout;
using std::endl;

typedef matrix<double, column_major, std::vector<double> > bmatrix;

void expm_najfeld_havel(bmatrix &AE);
void expm_higham05(bmatrix &AE);

int main () {

    bmatrix AE(4,4);

    AE(0,0) = 1.e+00; AE(0,1) = 1.e+02; AE(0,2) = 1.e+04; AE(0,3) = 1.e+06;
    AE(1,0) = 1.e-02; AE(1,1) = 1.e+00; AE(1,2) = 1.e+02; AE(1,3) = 1.e+04;
    AE(2,0) = 1.e-04; AE(2,1) = 1.e-02; AE(2,2) = 1.e+00; AE(2,3) = 1.e+02;
    AE(3,0) = 1.e-06; AE(3,1) = 1.e-04; AE(3,2) = 1.e-02; AE(3,3) = 1.e+00;

    cout << "\nA      = " << AE << endl;
    expm_najfeld_havel(AE);
    //expm_higham05(AE);
    cout << "exp(A) = " << AE << endl;

    AE(0,0) =-1.e+00; AE(0,1) =-2.e+00; AE(0,2) = 1.e+00; AE(0,3) = 0.e+00;
    AE(1,0) = 1.e+00; AE(1,1) = 1.e+00; AE(1,2) = 1.e+00; AE(1,3) = 0.e+00;
    AE(2,0) = 0.e+00; AE(2,1) = 0.e+00; AE(2,2) = 0.e+00; AE(2,3) = 0.e+00;
    AE(3,0) = 1.e+00; AE(3,1) = 1.e+00; AE(3,2) = 1.e+00; AE(3,3) = 0.e+01;

    cout << "\nA      = " << AE << endl;
    expm_najfeld_havel(AE);
    //expm_higham05(AE);
    cout << "exp(A) = " << AE << endl;

}

