#ifndef LINALG_H
#define LINALG_H

#include <string>

class Matrix; // forward decleration for decleration of friend function

class Vector {
private:
    double* elements;

public:
    int n;

    Vector(int, double);

    double& operator[](int);

    std::string to_string();

    double dot(Vector other);

    friend Vector multMatVec(double alpha, Matrix matA, Vector vecX, double beta, Vector vecY);
};

class Matrix {
private:
    double* elements;

public:
    int n, m;
    Matrix(int, int, double);

    double* operator[](int);

    std::string to_string();

    double* eigenValues();

    friend Vector multMatVec(double alpha, Matrix matA, Vector vecX, double beta, Vector vecY);
};

Vector multMatVec(double alpha, Matrix matA, Vector vecX, double beta, Vector vecY);

void TestMatrix();
void TestVector();
void TestFriendFunctions();

#endif
