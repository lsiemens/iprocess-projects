#ifndef GRAPH_SYSTEM_CPP
#define GRAPH_SYSTEM_CPP

#include "linalg.h"

#include <string>

void graph_spring_system(Vector vecX_x, Vector vecX_y, Matrix matK, std::string fname = "");
void graph_time_evolution(Vector vecT, Vector vecX);

#endif
