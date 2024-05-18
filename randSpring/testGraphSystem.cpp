#include "linalg.h"
#include "graph_system.h"
#include <iostream>

int main() {
    int numPoints = 10;
    Vector x = Vector(numPoints, 0);
    Vector y = Vector(numPoints, 10);
    Matrix adj = Matrix(numPoints, 0);

    double t = 0, dt = 1.0/numPoints;
    for (int i = 0; i < numPoints; i++) {
        t = i*dt;
        x[i] = t;
        y[i] = t*t;
    }
    adj[0][1] = 1;
    adj[1][0] = 1;

    adj[2][4] = 1;
    adj[4][2] = 1;

    adj[1][8] = 1;
    adj[8][1] = 1;

    adj[0][2] = 1;
    adj[2][0] = 1;

    graph(x, y, adj);
    return 0;
}
