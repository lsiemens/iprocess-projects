//
// Simple simulation of a mass on a spring. Print out the results for
// gnuplot to show the results.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void simulation(std::ofstream& fout, double m, double k, double x_0, double x_dot_0, double t_min, double t_max, int steps, int verbose=0) {
    // Initalize Arrays
    double dt = (t_max - t_min)/(steps - 1);

    std::vector<double> t(steps);
    for (int i = 0; i < steps; i++) {
        t[i] = t_min + dt*i;
    }

    std::vector<double> x(steps, x_0);
    std::vector<double> x_dot(steps, x_dot_0);
    std::vector<double> energy(steps, 0.0);
    energy[0] = 0.5*m*x_dot[0]*x_dot[0] + 0.5*k*x[0]*x[0];

    // Simulation

    for (int i = 1; i < steps; i++) {
        x_dot[i] = x_dot[i - 1] + dt*(-(k/m)*x[i - 1]);
        x[i] = x[i - 1] + dt*x_dot[i - 1];

        energy[i] = 0.5*m*x_dot[i]*x_dot[i] + 0.5*k*x[i]*x[i];

        if (verbose > 0) {
            std::cout << "x: " << x[i] << ", energy: " << energy[i] << std::endl;
        }
    }

    // Save to disk

    fout << "#t x x_dot energy" << std::endl;
    for (int i = 0; i < steps; i++) {
        fout << t[i] << " " << x[i] << " " << x_dot[i] << " " << energy[i] << std::endl;
    }
    fout << std::endl << std::endl;
    return;
}

int main() {
    std::ofstream fout("spring_simple.dat");
    simulation(fout, 1.0, 1.0, 1.0, 0.0, 0.0, 10.0, 1000);
    simulation(fout, 1.0, 1.0, 1.0, 0.0, 0.0, 100.0, 1000);
    fout.close();
}
