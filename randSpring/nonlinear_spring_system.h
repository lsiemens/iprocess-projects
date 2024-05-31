#ifndef NONLINEAR_SPRING_SYSTEM_H
#define NONLINEAR_SPRING_SYSTEM_H

#include "linalg.h"

#include <string>

class Nonlinear_Spring_System {
private:
    int num_particles;
    Vector vecM;
    Matrix matK;
    Matrix matL;
    Vector vecX_x0;
    Vector vecX_y0;
    Vector vecV_x0;
    Vector vecV_y0;

    Vector vecX_x;
    Vector vecX_y;
    Vector vecV_x;
    Vector vecV_y;
    double t;
    double energy;

    std::string tmp_dir = "tmp";

public:
    Nonlinear_Spring_System(int num_masses);

    void init_grid(double mass=1.0, double stiffness=1.0, double length=1.0);
    void init_gaussian(double mass=1.0, double stiffness=1.0, double length=1.0);
    void reset_simulation();

    void setup_tmp_dir();
    void animate_frame(std::string fname, int i, int i_max);

    void update_energy();
    void simulate_step(double dt);
    void simulation(std::string fname, double tMax, int steps);
};

#endif
