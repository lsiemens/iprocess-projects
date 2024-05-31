#include "graph_system.h"

#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include "gnuplot-iostream.h"

void graph_spring_system(Vector vecX_x, Vector vecX_y, Matrix matK, std::string fname) {
    Gnuplot gp;

    double x_0 = 0, y_0 = 0;
    for (int i = 0;  i < vecX_x.n; i++) {
        x_0 += vecX_x[i]/vecX_x.n;
        y_0 += vecX_y[i]/vecX_y.n;
    }

    double tmp_x = 0, tmp_y = 0;
    for (int i = 0;  i < vecX_x.n; i++) {
        tmp_x += (vecX_x[i] - x_0)*(vecX_x[i] - x_0);
        tmp_y += (vecX_y[i] - y_0)*(vecX_y[i] - y_0);
    }
    double x_std = std::sqrt(tmp_x/vecX_x.n);
    double y_std = std::sqrt(tmp_y/vecX_y.n);
    int width = 3;

    if (fname != "") {
//        gp << "set terminal png size 1920,1080\n";
        gp << "set terminal png size 840,680\n";
        gp << "set output '" + fname + "'\n";
    }
    gp << "set xrange [" << std::to_string(x_0 - width*x_std) << ":" << std::to_string(x_0 + width*x_std) << "]\n";
    gp << "set yrange [" << std::to_string(y_0 - width*y_std) << ":" << std::to_string(y_0 + width*y_std) << "]\n";
//    gp << "set yrange [-1:1]\n";
    gp << "plot " << gp.file1d(boost::make_tuple(vecX_x.to_vector(), vecX_y.to_vector())) << " with points pointtype 7 pointsize 2 notitle,\\\n";
    for (int i = 0; i < matK.n; i++) {
        for (int j = 0; j < i; j++) {
            if (matK[i][j] != 0.0) {
                std::vector<double> x = std::vector<double>(2);
                std::vector<double> y = std::vector<double>(2);
                x[0] = vecX_x[i];
                x[1] = vecX_x[j];
                y[0] = vecX_y[i];
                y[1] = vecX_y[j];
                gp << gp.file1d(boost::make_tuple(x, y)) << "with lines linecolor rgb 'black' linewidth 1 notitle,\\\n";
            }
        }
    }
    gp << "\n";
    gp << "pause mouse close\n";
}

void graph_time_evolution(Vector vecT, Vector vecX) {
    Gnuplot gp;
//    gp << "set xrange [-1:1]\n";
//    gp << "set yrange [-1:1]\n";
    gp << "plot " << gp.file1d(boost::make_tuple(vecT.to_vector(), vecX.to_vector())) << " with lines notitle\n";
    gp << "pause mouse close\n";
}
