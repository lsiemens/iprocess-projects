//
// Linear algebra for spring simulations
//
#include "linalg.h"

#include <iostream>
#include <string>

#include <cblas.h>
#include <lapacke.h>

// ------------ Vector Methods ------------
Vector::Vector() {
    n = 0;
    elements = {};
}

Vector::Vector(int size, double value=0) {
    n = size;
    elements = new double[size];
    for (int i = 0; i < n; i++) {
        elements[i] = value;
    }
}

double& Vector::operator[](int index) {
    return elements[index];
}

std::string Vector::to_string() {
    std::string result = "n=" + std::to_string(n) +  "\n[";
    for (int i = 0; i < n; i++) {
        result += std::to_string(elements[i]);
        if (i < n - 1) {
            result += ", ";
        }
    }
    result += "]";
    return result;
}

std::vector<double> Vector::to_vector() {
    return std::vector<double>(elements, elements + n);
}

double Vector::dot(Vector other) {
    return cblas_ddot(n, elements, 1, other.elements, 1);
}

// ------------ Matrix Methods ------------
Matrix::Matrix() {
    n = 0;
    m = 0;
    elements = {};
}

Matrix::Matrix(int size_n, double value = 0) {
    n = size_n;
    m = size_n;
    elements = new double[n*m];
    for (int i = 0; i < n*m; i++) {
        elements[i] = value;
    }
}

Matrix::Matrix(int size_n, int size_m, double value = 0) {
    n = size_n;
    m = size_m;
    elements = new double[n*m];
    for (int i = 0; i < n*m; i++) {
        elements[i] = value;
    }
}

double* Matrix::operator[](int index) {
    return (elements + m*index);
}

std::string Matrix::to_string() {
    std::string result = "n x m=" + std::to_string(n) + " x " + std::to_string(m) + "\n[[";
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            result += " [";
        }
        for (int j = 0; j < m; j++) {
            result += std::to_string(elements[m*i + j]);
            if (j < m - 1) {
                result += ", ";
            }
        }
        if (i < n - 1) {
            result += "],\n";
        }
    }
    result += "]]\n";
    return result;
}

double* Matrix::eigenValues() {
    double* eigenReal = new double[n];
    double* eigenImag = new double[n];
    int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, elements, n, eigenReal, eigenImag, NULL, 1, NULL, 1);

    std::cout << "DGEEV: info -> " << info << "\n";
    for (int i = 0; i < n; i++) {
        std::cout << "(" << eigenReal[i] << " + " << eigenImag[i] << "I), ";
    }
    return eigenReal;
}

// ------------ Friend Functions ------------
Vector multMatVec(double alpha, Matrix matA, Vector vecX, double beta, Vector vecY) {
    // Blas DGEMV solving y = alpha*A*x + beta*y
    cblas_dgemv(CblasColMajor, CblasNoTrans, matA.m, matA.n, alpha, matA.elements, matA.m, vecX.elements, 1, beta, vecY.elements, 1);
    return vecY;
}

// ---------------- TESTING ---------------- //
void TestMatrix() {
    Matrix matA = Matrix(3, 3, 0.0);
    for (int i = 0; i < matA.n; i++) {
        for (int j = 0; j < matA.m; j++) {
            matA[i][j] = (i + 1) + 100*(j + 1);
        }
    }

    for (int i = 0; i < matA.n; i++) {
        for (int j = 0; j < matA.m; j++) {
            std::cout << matA[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << matA.to_string() << "\n";

    Matrix matB = Matrix(2, 5, 0.0);
    for (int i = 0; i < matB.n; i++) {
        for (int j = 0; j < matB.m; j++) {
            matB[i][j] = (i + 1) + 100*(j + 1);
        }
    }

    for (int i = 0; i < matB.n; i++) {
        for (int j = 0; j < matB.m; j++) {
            std::cout << matB[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << matB.to_string() << "\n";
    Matrix matC = Matrix(3, 0, 0.0);
    std::cout << matC.to_string() << "\n";
    Matrix matD = Matrix(0, 3, 0.0);
    std::cout << matD.to_string() << "\n";
    Matrix matE = Matrix(0, 0, 0.0);
    std::cout << matE.to_string() << "\n";

    matA.eigenValues();
}

void TestVector() {
    // Test Vector element assignment and access
    Vector x = Vector(8);
    for (int i = 0; i < x.n; i++) {
        x[i] = i*i + 1;
    }

    for (int i = 0; i < x.n; i++) {
        std::cout << i << " " << x[i] << "\n";
    }

    std::cout << x.to_string() << "\n";
    std::cout << "x dot x = " << x.dot(x) << "\n";
    Vector y = Vector(0);
    std::cout << y.to_string() << "\n";
}

void TestFriendFunctions() {
    Matrix matA = Matrix(3, 3, 1.0);
    Vector x = Vector(3, 2.0);
    Vector y = Vector(3, 3.0);
    std::cout << "\nmatA: " << matA.to_string() << "\n";
    std::cout << "x: " << x.to_string() << "\n\n";
    std::cout << "y: " << y.to_string() << "\n\n";

    Vector z = multMatVec(3.0, matA, x, 4.0, y);
    std::cout << "z: " << z.to_string() << "\n\n";
}
