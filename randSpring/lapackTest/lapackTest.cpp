//
// Test finding eigenvalues with LAPACK
//

#include <iostream>
#include <fstream>

using namespace std;

// Use dgeev_ from LAPACK
extern "C" {
    // dgeev_ stands for Double GEneral matrix EigenValue. It solves for
    // the eigen values of an arbitrary matrix of doubles.

    // char (in): JOBVL
    // char (in): JOBVR

    // int            (in): N
    // double[][] (in/out): A
    // int            (in): LDA

    // double[] (out): WR
    // double[] (out): WI

    // double[][] (out): VL
    // int         (in): LDVL
    // double[][] (out): VR
    // int         (in): LDVR

    // double[] (out): WORK
    // int       (in): LWORK

    // int (out): INFO
    extern int dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

int main() {
    int n;
    cout << "Enter matrix size: ";
    cin >> n;

    double data[n*n];

    cout << "Enter matrix values!" << endl;
    for (int i=0; i < n; i++) {
        cout << "Row: " << i << endl;
        for (int j=0; j < n; j++) {
            cout << "Column: " << j << endl;
            cout << "Enter matrix element: ";
            cin >> data[i + n*j];
        }
    }

    cout << endl << "Your matrix is," << endl << "[";
    for (int i=0; i < n; i++) {
        if (i > 0) {
            cout << " ";
        }
        cout << "[ ";
        for (int j=0; j < n; j++) {
            cout << data[i + n*j] << " ";
        }
        cout << "]";
        if (i < n - 1) {
            cout << endl;
        }
    }
    cout << "]" << endl;

    char Nchar='N';
    double* eigReal = new double[n];
    double* eigImag = new double[n];
    double* vl;
    double* vr;
    int one=1;
    int lwork=3*n;
    double* work=new double[lwork];
    int info;

    dgeev_(&Nchar, &Nchar, &n, data, &n, eigReal, eigImag, vl, &one, vr, &one, work, &lwork, &info);

    if (info != 0) {
        cout << "DGEEV: Failed with error code -> " << info << endl;
    }

    cout << "EIGENVALUES" << endl;
    for (int i=0; i<n; i++) {
        cout << eigReal[i] << "+" << eigImag[i] << "J, ";
    }
    cout <<  endl;
}
