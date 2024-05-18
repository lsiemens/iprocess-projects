//
// Simulation of a 2D system of springs. Use the gnuplot library to show
// the results.
//
// M is a vector of masses m_i
// K is a symmetric matrix of spring constants with K_{ij}=0 along the diagonal
// L is a symmetric matrix of spring lengths with L_{ij}=0 along the diagonal
//
// In 2D the "vector" x_k is a vector of 2D vectors. These will be represented
// as two vectors x_{xi} and x_{yi}, Ill use the index "a" for these components.
//
// y is the two matricies matrix y_{aij} = x_{ai} - x_{aj} - L_{ij}(x_{ai} - x_{aj})/|x_i - x_j|
// Note, the matrices y are anitisymmetric and nonlinear functions of x_{ai}.
//
// The equation of motion is a_{ak} = (1/M_k) SUM_i K_{ik}y_{aik}, where a is x_dot_dot
//
// Simulate using the semi-implicit Euler integrator:
// v_{i+1} = v_i + dt*A*x_i
// x_{i+1} = x_i + dt*v_i
//

#include "nonlinear_spring_system.h"

#include <iostream>
#include <cmath>

#include "graph_system.h"

Nonlinear_Spring_System::Nonlinear_Spring_System(int grid_size, double mass, double stiffness, double length) {
    num_particles = grid_size*grid_size;
    vecM = Vector(num_particles, mass);
    matK = Matrix(num_particles, 0);
    matL = Matrix(num_particles, 0);
    vecX_x0 = Vector(num_particles, 0.0);
    vecX_y0 = Vector(num_particles, 0.0);
    vecV_x0 = Vector(num_particles, 0.0);
    vecV_y0 = Vector(num_particles, 0.0);

    int i = 0, j = 0;
    for (int grid_i = 0; grid_i < grid_size; grid_i++) {
        for (int grid_j = 0; grid_j < grid_size; grid_j++) {
            i = grid_i + grid_j*grid_size;
            vecX_x0[i] = grid_i*length;
            vecX_y0[i] = grid_j*length;

            if (grid_i > 0) {
                j = (grid_i - 1) + grid_j*grid_size;
                matK[i][j] = stiffness;
                matK[j][i] = stiffness;
                matL[i][j] = length;
                matL[j][i] = length;
            }

            if (grid_j > 0) {
                j = grid_i + (grid_j - 1)*grid_size;
                matK[i][j] = stiffness;
                matK[j][i] = stiffness;
                matL[i][j] = length;
                matL[j][i] = length;
            }
        }
    }

//    vecM[0] = INFINITY; // make infine mass so the spring system does not drift away

    std::cout << matK.to_string() << "\n";
    std::cout << vecX_x0.to_string() << "\n";
    std::cout << vecX_y0.to_string() << "\n";
    std::cout << vecM.to_string() << "\n";
    graph(vecX_x0, vecX_y0, matK);
}

int main() {
    Nonlinear_Spring_System sys = Nonlinear_Spring_System(25);
    return 0;
}
