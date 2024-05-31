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
#include <filesystem>
#include <cmath>
#include <random>
#include <ctime>

#include "graph_system.h"

Nonlinear_Spring_System::Nonlinear_Spring_System(int num_masses) {
    num_particles = num_masses;
    vecM = Vector(num_particles, INFINITY);
    matK = Matrix(num_particles, 0.0);
    matL = Matrix(num_particles, 0.0);
    vecX_x0 = Vector(num_particles, 0.0);
    vecX_y0 = Vector(num_particles, 0.0);
    vecV_x0 = Vector(num_particles, 0.0);
    vecV_y0 = Vector(num_particles, 0.0);
}

void Nonlinear_Spring_System::init_grid(double mass, double stiffness, double length) {
    int i = 0, j = 0;
    int grid_size = static_cast<int>(std::floor(std::sqrt(num_particles)));
    for (int grid_i = 0; grid_i < grid_size; grid_i++) {
        for (int grid_j = 0; grid_j < grid_size; grid_j++) {
            i = grid_i + grid_j*grid_size;
            vecM[i] = mass;
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

            if ((grid_i > 0) && (grid_j > 0)) {
                j = (grid_i - 1) + (grid_j - 1)*grid_size;
                matK[i][j] = stiffness;
                matK[j][i] = stiffness;
                matL[i][j] = length;
                matL[j][i] = length;
            }

            if ((grid_i > 0) && (grid_j + 1 < grid_size)) {
                j = (grid_i - 1) + (grid_j + 1)*grid_size;
                matK[i][j] = stiffness;
                matK[j][i] = stiffness;
                matL[i][j] = length;
                matL[j][i] = length;
            }
        }
    }

    //vecM[0] = INFINITY; // make infine mass so the spring system does not drift away
    //vecM[num_particles - 1] = INFINITY; // make infine mass so the spring system does not drift away

    vecX_x0[2] = 1.8;
    vecX_y0[2] = 0.2;
    //vecV_x0[num_particles - 2] = 0.1;
    //vecV_y0[num_particles - 2] = 0.1;
    //vecV_x0[1] = -0.1;
    //vecV_y0[1] = -0.1;

    std::cout << matK.to_string() << "\n";
    std::cout << vecX_x0.to_string() << "\n";
    std::cout << vecX_y0.to_string() << "\n";
    std::cout << vecM.to_string() << "\n";
    graph_spring_system(vecX_x0, vecX_y0, matK);

    reset_simulation();
}

void Nonlinear_Spring_System::init_gaussian(double mass, double stiffness, double length) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> standard_normal(0.0, 1.0);
    std::uniform_real_distribution<double> standard_uniform(0.0, 1.0);

    for (int i = 0; i < num_particles; i++) {
        vecM[i] = mass;
        vecX_x0[i] = standard_normal(generator);
        vecX_y0[i] = standard_normal(generator);

        for (int j = 0; j < i; j++) {
            if (standard_uniform(generator) < 0.6) {
                matK[i][j] = std::abs(standard_normal(generator));
                matK[j][i] = matK[i][j];
                matL[i][j] = 1.0;
                matL[j][i] = matL[i][j];

//                matK[i][j] = std::pow(standard_normal(generator), 2);
//                matK[j][i] = matK[i][j];
//                matL[i][j] = std::pow(standard_normal(generator) + 0.5, 2);
//                matL[j][i] = matL[i][j];
            }
        }
    }

    vecM[0] = INFINITY; // make infine mass so the spring system does not drift away
    //vecM[num_particles - 1] = INFINITY; // make infine mass so the spring system does not drift away

    //vecX_x0[2] = 1.8;
    //vecX_y0[2] = 0.2;
    //vecV_x0[num_particles - 2] = 0.1;
    //vecV_y0[num_particles - 2] = 0.1;
    //vecV_x0[1] = -0.1;
    //vecV_y0[1] = -0.1;

    std::cout << matK.to_string() << "\n";
    std::cout << vecX_x0.to_string() << "\n";
    std::cout << vecX_y0.to_string() << "\n";
    std::cout << vecM.to_string() << "\n";
    graph_spring_system(vecX_x0, vecX_y0, matK);

    reset_simulation();
}

void Nonlinear_Spring_System::reset_simulation() {
    vecX_x = Vector(num_particles, 0);
    vecX_y = Vector(num_particles, 0);
    vecV_x = Vector(num_particles, 0);
    vecV_y = Vector(num_particles, 0);
    for (int i = 0; i < num_particles; i++) {
        vecX_x[i] = vecX_x0[i];
        vecX_y[i] = vecX_y0[i];

        vecV_x[i] = vecV_x0[i];
        vecV_y[i] = vecV_y0[i];
    }
    t = 0;
    update_energy();
}

void Nonlinear_Spring_System::setup_tmp_dir() {
    if (!std::filesystem::exists(tmp_dir)) {
        std::filesystem::create_directory(tmp_dir);
    }
}

void Nonlinear_Spring_System::animate_frame(std::string fname, int i, int i_max) {
    // Fixed width for the index
    int width = std::to_string(i_max).size();
    std::string index = std::to_string(i);
    std::string frame_name = fname + "_" + std::string(width - index.size(), '0') + index + ".png";
    graph_spring_system(vecX_x, vecX_y, matK, tmp_dir + "/" + frame_name);
}

void Nonlinear_Spring_System::update_energy() {
    energy = 0;
    for (int i = 0; i < num_particles; i++) {
        if (vecM[i] != INFINITY) {
            // Add kinetic energy for particle i if it is not a fixed point
            energy += 0.5*vecM[i]*(vecV_x[i]*vecV_x[i] + vecV_y[i]*vecV_y[i]);
        }

        for (int j = 0; j < i; j++) {
            double xx_xij = vecX_x[i] - vecX_x[j];
            double xx_yij = vecX_y[i] - vecX_y[j];
            double normxx_ij = std::sqrt(xx_xij*xx_xij + xx_yij*xx_yij);
            // Add potential energy for spring i,j
            energy += 0.5*matK[i][j]*(normxx_ij - matL[i][j])*(normxx_ij - matL[i][j]);
        }
    }
}

void Nonlinear_Spring_System::simulate_step(double dt) {
    for (int k = 0; k < num_particles; k++) {
        double SUM_K_Y_xk = 0;
        double SUM_K_Y_yk = 0;
        for (int i = 0; i < num_particles; i++) {
            double xx_xik = vecX_x[i] - vecX_x[k];
            double xx_yik = vecX_y[i] - vecX_y[k];
            // Add very small positive number to the norm to regularize the value
            double normxx_ik = std::sqrt(xx_xik*xx_xik + xx_yik*xx_yik + 1E-32);
            SUM_K_Y_xk += matK[i][k]*(xx_xik - matL[i][k]*xx_xik/normxx_ik);
            SUM_K_Y_yk += matK[i][k]*(xx_yik - matL[i][k]*xx_yik/normxx_ik);
        }
        // forward Euler
        //double tmpV_x = vecV_x[k];
        //double tmpV_y = vecV_y[k];
        vecV_x[k] = vecV_x[k] + dt*SUM_K_Y_xk/vecM[k] - dt*0.8*vecV_x[k];
        vecV_y[k] = vecV_y[k] + dt*SUM_K_Y_yk/vecM[k] - dt*0.8*vecV_y[k];

        //vecX_x[k] = vecX_x[k] + dt*tmpV_x;
        //vecX_y[k] = vecX_y[k] + dt*tmpV_y;
        vecX_x[k] = vecX_x[k] + dt*vecV_x[k];
        vecX_y[k] = vecX_y[k] + dt*vecV_y[k];
    }

    t += dt;
    update_energy();
}

void Nonlinear_Spring_System::simulation(std::string fname, double tMax, int steps) {
    setup_tmp_dir();
    double dt = tMax/(steps - 1);

    Vector vecT = Vector(steps, 0);
    Vector vecE = Vector(steps, 0);
    for (int i = 0; i < steps; i++) {
        if (i == 0) {
            // Initial configuration
        } else {
            simulate_step(dt);
        }
        std::cout << "\rSimulation | Step: [" << i + 1 << "/" << steps << "], ";
        std::cout << "Time: [" << t << "/" << tMax << "], ";
        std::cout << "Energy: " << energy;
        std::cout.flush();
        vecT[i] = t;
        vecE[i] = energy;
        animate_frame(fname, i, steps);
    }
    std::cout << "\n";
    graph_time_evolution(vecT, vecE);
    graph_spring_system(vecX_x, vecX_y, matK);
}

int main() {
    Nonlinear_Spring_System sys = Nonlinear_Spring_System(20);
    sys.init_gaussian();
    sys.simulation("test", 10.0, 1000);
    return 0;
}
