//
// Simulation of a mass on a spring. Use the gnuplot library to show
// the results.
//
// Simulate using the semi-implicit Euler integrator:
// x_dot_{i+1} = x_dot_i + dt*(-k)*x_i
// x_{i+1} = x_i + dt*x_dot_{i+1}/m
//

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"

void simulation(Gnuplot& gp, double m, double k, double x_0, double x_dot_0, double t_min, double t_max, int steps, int verbose=0) {
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
        x_dot[i] = x_dot[i - 1] + dt*(-k/m)*x[i - 1];
        x[i] = x[i - 1] + dt*x_dot[i];

        energy[i] = 0.5*m*x_dot[i]*x_dot[i] + 0.5*k*x[i]*x[i];

        if (verbose > 0) {
            std::cout << "x: " << x[i] << ", energy: " << energy[i] << std::endl;
        }
    }

    // Add plot with Gnuplot

    gp << "plot" << gp.file1d(boost::make_tuple(t, x)) << "with lines title 'x',\\" << std::endl;
    gp << gp.file1d(boost::make_tuple(t, energy)) << "with lines title 'energy'" << std::endl;
    return;
}

int main() {
    // Graph with Gnuplot

    Gnuplot gp;
    gp << "set xlabel 'Time t'" << std::endl;
    gp << "set multiplot layout 1,2 columns title 'Symplectic Spring Simulation'" << std::endl;
    gp << "set title 'Short duration'" << std::endl;
    simulation(gp, 1.0, 1.0, 1.0, 0.0, 0.0, 10.0, 1000);
    gp << "set title 'Long duration'" << std::endl;
    simulation(gp, 1.0, 1.0, 1.0, 0.0, 0.0, 100.0, 1000);
    gp << "unset multiplot" << std::endl;
    return 0;
}
