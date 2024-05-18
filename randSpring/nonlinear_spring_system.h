#ifndef NONLINEAR_SPRING_SYSTEM_H
#define NONLINEAR_SPRING_SYSTEM_H

#include "linalg.h"

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

public:
    Nonlinear_Spring_System(int grid_size, double mass=1.0, double stiffness=1.0, double length=1.0);

    void simulation(double tMax, int steps);
};

#endif
