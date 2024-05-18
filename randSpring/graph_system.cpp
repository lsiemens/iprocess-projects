#include "graph_system.h"

#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include "gnuplot-iostream.h"

void graph(Vector vecX_x, Vector vecX_y, Matrix matK) {
    Gnuplot gp;
//    gp << "set xrange [-1:1]\n";
//    gp << "set yrange [-1:1]\n";
    gp << "plot " << gp.file1d(boost::make_tuple(vecX_x.to_vector(), vecX_y.to_vector())) << " with points pointtype 7 pointsize 3 notitle,\\\n";
    for (int i = 0; i < matK.n; i++) {
        for (int j = 0; j < i; j++) {
            if (matK[i][j] != 0.0) {
                std::vector<double> x = std::vector<double>(2);
                std::vector<double> y = std::vector<double>(2);
                x[0] = vecX_x[i];
                x[1] = vecX_x[j];
                y[0] = vecX_y[i];
                y[1] = vecX_y[j];
                gp << gp.file1d(boost::make_tuple(x, y)) << "with lines linecolor rgb 'black' linewidth 4 dashtype '--' notitle,\\\n";
            }
        }
    }
    gp << "\n";
    gp << "pause mouse close\n";
}
