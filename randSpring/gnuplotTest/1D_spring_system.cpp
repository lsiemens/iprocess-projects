//
// Simulation of a 1D system of springs. Use the gnuplot library to show
// the results.
//
// M is a diagonal matrix of the masses
// K is a symmetric matrix of spring constants with K_{ij}=0 along the diagonal
//
// K' is the matrix with elements such that, if (i!=j) then K'_{ij} = -K_{ij}
// and if (i==j) then K'_{ii} = sum_k K_{ik}
//
// The equation of motion is Ma = -K'x, where a is x_dot_dot
//
// In this 1d system all masses have mass m, except for the first and
// last which are infinite. All of the masses are connected to their
// neighbor on each side with a spring with spring constant k. so
//
//     ┌─         ─┐
//     │ ∞ 0 0 . . │
//     │ 0 m 0     │
// M = │ 0 0 m     │
//     │ .     .   │
//     │ .       . │
//     └─         ─┘
//
//     ┌─         ─┐
//     │ 0 k 0 0 . │
//     │ k 0 k 0   │
// K = │ 0 k 0 k   │
//     │ 0 0 k 0   │
//     │ .       . │
//     └─         ─┘
//
// Then the matricies M^-1 and K' are,
//
//        ┌─                   ─┐
//        │  0   0   0   .   .  │
//        │  0  1/m  0          │
// M^-1 = │  0   0  1/m         │
//        │  .           .      │
//        │  .               .  │
//        └─                   ─┘
//
//      ┌─              ─┐
//      │  k -k  0  0  . │
//      │ -k 2k -k  0    │
// K' = │  0 -k 2k -k    │
//      │  0  0 -k 2k    │
//      │  .           . │
//      └─              ─┘
//
// Finally define A = -M^-1 * K'
//
//     ┌─                       ─┐
//     │   0     0     0     .   │
//     │   k/m -2k/m   k/m       │
// A = │   0     k/m -2k/m       │
//     │   0     0     k/m       │
//     │   .                 .   │
//     └─                       ─┘
//
// Simulate using the semi-implicit Euler integrator:
// x_dot_{i+1} = x_dot_i + dt*A*x_i
// x_{i+1} = x_i + dt*x_dot_i
//

#include <stdexcept>
#include <filesystem>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"

// Function definitions
void printMatrix(std::string msg, std::vector<double> vectorV);
void printMatrix(std::string msg, std::vector<std::vector<double>> matrixM);
std::vector<double> multMatrix(std::vector<std::vector<double>> matrixM, std::vector<double> vectorV);
double multMatrix(std::vector<double> vectorV, std::vector<double> vectorW);
std::vector<double> multMatrix(std::vector<double> vectorV, double scalarS);
std::vector<double> addMatrix(std::vector<double> vectorV, std::vector<double> vectorW);
void simulation(std::string fname, double m, double k, std::vector<double>& x_0, std::vector<double>& x_dot_0, double t_min, double t_max, int steps, int verbose=0);
void animation(std::string fname, std::vector<double> t, std::vector<std::vector<double>> x, std::vector<double> energy);
void renderFrame(std::string fname, double t, std::vector<double> x, double energy);

// Linear algebra functions and utilities
void printMatrix(std::string msg, std::vector<double> vectorV) {
    // Print Vector
    std::cout << msg << " Size=" << vectorV.size() << ", [ ";
    for (double v: vectorV) {
        std::cout << v << " ";
    }
    std::cout << "]" << std::endl;
    return;
}

void printMatrix(std::string msg, std::vector<std::vector<double>> matrixM) {
    // Print Matrix
    std::cout << msg << " Size=" << matrixM.size() << "x" << matrixM[0].size() << ", [ ";
    for (std::vector<double> row: matrixM) {
        std::cout << "[";
        for (double m: row) {
            std::cout << m << " ";
        }
        std::cout << "] ";
    }
    std::cout << "]" << std::endl;
    return;
}

std::vector<double> multMatrix(std::vector<std::vector<double>> matrixM, std::vector<double> vectorV) {
    // Matrix Vector Multiplication
    int sizeMn = matrixM.size(), sizeMm = matrixM[0].size(), sizeVn = vectorV.size();
    // Sanity check
    if (sizeMm != sizeVn) {
        throw std::invalid_argument(std::to_string(sizeMn) + "x" + std::to_string(sizeMm) + " matrix can not be multiplied with a " + std::to_string(sizeVn) + " dimentional vector!");
    }

    std::vector<double> result(sizeMn);
    for (int i=0; i < sizeMn; i++) {
        result[i] = multMatrix(matrixM[i], vectorV);
    }
    return result;
}

double multMatrix(std::vector<double> vectorV, std::vector<double> vectorW) {
    // Vector, Dot Product
    int sizeVn = vectorV.size(), sizeWn = vectorW.size();
    // Sanity check
    if (sizeVn != sizeWn) {
        throw std::invalid_argument("Can not take the dot product of a " + std::to_string(sizeVn) + " dimentional vector with a " + std::to_string(sizeWn) + " dimentional vector!");
    }

    double result = 0;
    for (int i=0; i < sizeVn; i++) {
        result += vectorV[i]*vectorW[i];
    }
    return result;
}

std::vector<double> multMatrix(std::vector<double> vectorV, double scalarS) {
    // Vector, Scalar Multiplication
    int sizeVn = vectorV.size();
    std::vector<double> result(sizeVn);
    for (int i=0; i < sizeVn; i++) {
        result[i] = scalarS*vectorV[i];
    }
    return result;
}

std::vector<double> addMatrix(std::vector<double> vectorV, std::vector<double> vectorW) {
    // Vector Addition
    int sizeVn = vectorV.size(), sizeWn = vectorW.size();
    // Sanity check
    if (sizeVn != sizeWn) {
        throw std::invalid_argument("Can not add a " + std::to_string(sizeVn) + " dimentional vector with a " + std::to_string(sizeWn) + " dimentional vector!");
    }

    std::vector<double> result(sizeVn);
    for (int i=0; i < sizeVn; i++) {
        result[i] = vectorV[i] + vectorW[i];
    }
    return result;
}

// Simulation
void simulation(std::string fname, double m, double k, std::vector<double>& x_0, std::vector<double>& x_dot_0, double t_min, double t_max, int steps, int verbose) {
    // Sanity check
    int num_mass = x_0.size();
    std::cout << "Number of particles: " << num_mass << std::endl;
    if (x_0.size() != x_dot_0.size()) {
        throw std::invalid_argument("Vectors x_0 and x_dot_0 must have the same length!");
    }

    // Initalize Dynamic Arrays
    double dt = (t_max - t_min)/(steps - 1);

    std::vector<double> t(steps);
    for (int i = 0; i < steps; i++) {
        t[i] = t_min + dt*i;
    }

    std::vector<std::vector<double>> x(steps, x_0);
    std::vector<std::vector<double>> x_dot(steps, x_dot_0);
    std::vector<double> energy(steps, 0.0);
    energy[0] = 0.5*m*multMatrix(x_dot[0], x_dot[0]) + 0.5*k*multMatrix(x[0], x[0]);

    // Initalize System Configuration Array
    std::vector<std::vector<double>> matrixA(num_mass, std::vector<double>(num_mass));
    for (int i = 1; i < num_mass - 1; i++) {
        matrixA[i][i - 1] = k/m;
        matrixA[i][i] = -2*k/m;
        matrixA[i][i + 1] = k/m;
    }

    // DEBUG
    printMatrix("A:", matrixA);

    // Simulation
    for (int i = 1; i < steps; i++) {
        x_dot[i] = addMatrix(x_dot[i - 1], multMatrix(multMatrix(matrixA, x[i - 1]), dt));
        x[i] = addMatrix(x[i - 1], multMatrix(x_dot[i], dt));

        energy[i] = 0.5*m*multMatrix(x_dot[i], x_dot[i]) + 0.5*k*multMatrix(x[i], x[i]);

        if (verbose > 0) {
            std::cout << "Step: " << i << std::endl << "t_i: " << t[i] << std::endl;
            printMatrix("x_i:", x[i]);
            printMatrix("x_dot_i:", x_dot[i]);
            std::cout << "energy_i: " << energy[i] << std::endl << std::endl;
        }
    }

    // Graph the energy using Gnuplot
    Gnuplot gp;
    gp << "set xlabel 'Time t'" << std::endl;
    gp << "set ylabel 'Energy'" << std::endl;
    gp << "set title 'Energy Evolution'" << std::endl;
    gp << "plot" << gp.file1d(boost::make_tuple(t, energy)) << "with lines title 'energy'" << std::endl;

    animation(fname, t, x, energy);
    return;
}

void animation(std::string fname, std::vector<double> t, std::vector<std::vector<double>> x, std::vector<double> energy) {
    using namespace std::filesystem;

    std::string tmp_dir = "tmp";

    // Creat temp directory if it does not exist
    if (!exists(tmp_dir)) {
        create_directory(tmp_dir);
    }

    int width = std::to_string(t.size() - 1).size(); // get fix with for index
    for (int i = 0; i < t.size(); i++) {
        std::string frame_name = std::to_string(i);
        frame_name = fname + "_" + std::string(width - frame_name.size(), '0') + frame_name + ".png";

        renderFrame(tmp_dir + "/" + frame_name, t[i], x[i], energy[i]);
        std::cout << "\rRender Frames: " << i << "/" << t.size();
        std::cout.flush();
    }
}

void renderFrame(std::string fname, double t, std::vector<double> x, double energy) {
    // Graph the springs using Gnuplot
    Gnuplot gp;
    gp << "set terminal png" << std::endl;
    gp << "set output '" + fname + "'" << std::endl;
    gp << "set title 'Spring Evolution'" << std::endl;
    gp << "set label 1 'Time t:" + std::to_string(t) + " Energy: " + std::to_string(energy) + "' at 0.5,0.5 offset 0,0" << std::endl;
    gp << "set yrange [-1:1]" << std::endl;
    gp << "set xrange [" << x[0] << ":" << x[x.size() - 1] <<"]" << std::endl;
    gp << "set xlabel 'Position x'" << std::endl;
    gp << "plot" << gp.file1d(boost::make_tuple(x, std::vector<double>(x.size(), 0))) << " with points pointtype 5 title 'Spring System',\\" << std::endl;
    gp <<           gp.file1d(boost::make_tuple(x, std::vector<double>(x.size(), 0))) << " with lines dashtype '--' notitle" << std::endl;
    return;
}

int main() {
    // Graph with Gnuplot
    int num_mass = 25;

    // Initalize inital values
    std::vector<double> x_0(num_mass);
    for (int i = 0; i < num_mass; i++) {
        x_0[i] = 1.0*i;
    }
    std::vector<double> x_dot_0(num_mass, 0.0);

    x_dot_0[1] = 3.0;

    // DEBUG
    printMatrix("x_0:", x_0);
    printMatrix("x_dot_0:", x_dot_0);

    simulation("1D_spring_system", 1.0, 1.0, x_0, x_dot_0, 0.0, 100.0, 1000);
    return 0;
}
