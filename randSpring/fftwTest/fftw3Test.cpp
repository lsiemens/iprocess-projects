//
// Test finding the Fourier transform of data
//

#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>

#include <complex>
#include <fftw3.h>

#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"

void graph(std::vector<double> x, std::vector<double> realy, std::vector<double> imagy, std::string xlabel) {
    std::vector<double> magy(x.size(), 0);
    for (int i = 0; i < x.size(); i++) {
        magy[i] = std::sqrt(realy[i]*realy[i] + imagy[i]*imagy[i]);
    }

    Gnuplot gp;
    gp << "set xlabel '" << xlabel << ",\n";
    gp << "plot" << gp.file1d(boost::make_tuple(x, realy)) << "with lines title 'Real',\\\n";
    gp << gp.file1d(boost::make_tuple(x, imagy)) << "with lines title 'Imag',\\\n";
    gp << gp.file1d(boost::make_tuple(x, magy)) << "with lines title 'Mag'\n";
    gp << "pause mouse close\n";
}

int main() {
    int sizeN = 40000;

    std::vector<double> t(sizeN, 0);
    std::vector<double> f(sizeN, 0);
    std::vector<double> realDataIn(sizeN, 0);
    std::vector<double> imagDataIn(sizeN, 0);
    std::vector<double> realDataOut(sizeN, 0);
    std::vector<double> imagDataOut(sizeN, 0);

    // Initialize t
    double tMax = 200.0;
    double dt = tMax/(sizeN - 1);
    for (int i = 0; i < sizeN; i++) {
        t[i] = -tMax/2 + dt*i;
    }

    // Initialize f
    for (int i = 0; i < sizeN; i++) {
        f[i] = -1/(2*dt) + i/(tMax);
    }

    // Initialize data
    for (int i = 0; i < sizeN; i++) {
        realDataIn[i] = std::cos(M_PI*M_PI*t[i]);
        realDataIn[i] += std::cos(2*10*M_PI*t[i]);
        realDataIn[i] += std::cos(2*M_E*M_PI*t[i]);
        realDataIn[i] *= std::pow(M_E, -std::pow(t[i]/7, 4)); // window
    }

    // Show Input data
    graph(t, realDataIn, imagDataIn, "t");

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sizeN);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sizeN);

    p = fftw_plan_dft_1d(sizeN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//    p = fftw_plan_dft_1d(sizeN, in, out, FFTW_FORWARD, FFTW_MEASURE);

    for (int i = 0; i < sizeN; i++) {
        int sign = 2*(i % 2) - 1;
        in[i][0] = sign*realDataIn[i];
        in[i][1] = sign*imagDataIn[i];
    }

    fftw_execute(p); // Transform what ever is in "in"

    for (int i = 0; i < sizeN; i++) {
        int sign = 2*(i % 2) - 1;
        realDataOut[i] = sign*out[i][0];
        imagDataOut[i] = sign*out[i][1];
    }

    // Show Input data
    graph(f, realDataOut, imagDataOut, "f");

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}
