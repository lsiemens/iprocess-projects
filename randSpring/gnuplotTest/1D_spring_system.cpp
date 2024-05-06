//
// Simulation of a 1D system of springs. Use the gnuplot library to show
// the results.
//
// M is a diagonal matrix of the masses
// K is a symmetric matrix of spring constants with K_{ij}=0 along the diagonal
//
// K' is the matrix with elements such that, if (i!=j) then K'_{ij} = -K_{ij}
// and if (i==j) then K'_{ii} = sum_k K_{ik}
//
// The equation of motion is Ma = -K'x, where a is x_dot_dot
//
// In this 1d system all masses have mass m, except for the first and
// last which are infinite. All of the masses are connected to their
// neighbour on each side with a spring with spring constant k. so
//
//     ┌─         ─┐
//     │ ∞ 0 0 . . │
//     │ 0 m 0     │
// M = │ 0 0 m     │
//     │ .     .   │
//     │ .       . │
//     └─         ─┘
//
//     ┌─         ─┐
//     │ 0 k 0 0 . │
//     │ k 0 k 0   │
// K = │ 0 k 0 k   │
//     │ 0 0 k 0   │
//     │ .       . │
//     └─         ─┘
//
// Then the matricies M^-1 and K' are,
//
//        ┌─                   ─┐
//        │  0   0   0   .   .  │
//        │  0  1/m  0          │
// M^-1 = │  0   0  1/m         │
//        │  .           .      │
//        │  .               .  │
//        └─                   ─┘
//
//      ┌─              ─┐
//      │  k -k  0  0  . │
//      │ -k 2k -k  0    │
// K' = │  0 -k 2k -k    │
//      │  0  0 -k 2k    │
//      │  .           . │
//      └─              ─┘
//
// Finally define A = -M^-1 * K'
//
//     ┌─                       ─┐
//     │   0     0     0     .   │
//     │   k/m -2k/m   k/m       │
// A = │   0     k/m -2k/m       │
//     │   0     0     k/m       │
//     │   .                 .   │
//     └─                       ─┘
//
// Simulate using the semi-implicit Euler integrator:
// x_dot_{i+1} = x_dot_i + dt*A*x_i
// x_{i+1} = x_i + dt*x_dot_i
//
// The distribution of frequencies in the system is found by analyzing the
// matricies. Using LAPACK to calculate the eigenvalues of the matrix
// A, ie -M^-1*K, the eigenvalues are all negative real numbers. This
// implies that solutions to the equations of motion are sinusoidal. The
// possible frequencies are given by the square root of the absolute value
// of the eigenvalues of the matrix A. Making a histogram of these
// frequencies then gives the distribution of frequencies of this system.
//
// the power spectrum for the system is found by using the Fourier transform
// of the x(t) of the masses. Using FFTW3 the power spectral density for
// the position of each (finite mass) particle is calculated, the average
// spectrum is taken to be the spectrum for the system. If the system is
// perturbed from equilibrium by kicking a single mass and allowed to evolve
// for much longer than the signal crossing time of the system, then the
// power spectrum a general trend of decreasing power for increasing
// frequencies. Note, that the density of spectral lines (if resolved)
// is increasing with increasing frequency up until a cutoff frequency
// F. This trend of the density of allowed frequencies, and maximum
// frequency F are both features shared with the histogram of the frequencies
// found by analyzing the eigenvalues of the system.
//

#include <stdexcept>
#include <filesystem>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include <complex>
#include <fftw3.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include "gnuplot-iostream.h"

// Function definitions
void printMatrix(std::string msg, std::vector<double> vectorV);
void printMatrix(std::string msg, std::vector<std::vector<double>> matrixM);
std::vector<double> multMatrix(std::vector<std::vector<double>> matrixM, std::vector<double> vectorV);
double multMatrix(std::vector<double> vectorV, std::vector<double> vectorW);
std::vector<double> multMatrix(std::vector<double> vectorV, double scalarS);
std::vector<double> addMatrix(std::vector<double> vectorV, std::vector<double> vectorW);

void eigenvalues(std::vector<std::vector<double>> matrixA);
void powerSpectrum(std::vector<double> t, std::vector<std::vector<double>> data);

void showPlot(std::vector<double> t, std::vector<double> y, std::string title="", std::string xlabel="", std::string ylabel="", std::string graphlabel="", double xMin=0, double xMax=0, double yMin=0, double yMax=0);
void plotHistogram(std::vector<double> vectorV, int bins);
void animation(std::string fname, std::vector<double> t, std::vector<std::vector<double>> x, std::vector<double> energy, int stride=1);
void renderFrame(std::string fname, double t, std::vector<double> x, double energy);

void simulation(std::string fname, double m, double k, std::vector<double>& x_0, std::vector<double>& x_dot_0, double tMax, int steps, int verbose=0);

// Use dgeev_ from LAPACK
extern "C" {
    // dgeev_ stands for Double GEneral matrix EigenValue. It solves for
    // the eigenvalues of an arbitrary matrix of doubles.

    // char (in): JOBVL
    // char (in): JOBVR

    // int            (in): N
    // double[][] (in/out): A
    // int            (in): LDA

    // double[] (out): WR
    // double[] (out): WI

    // double[][] (out): VL
    // int         (in): LDVL
    // double[][] (out): VR
    // int         (in): LDVR

    // double[] (out): WORK
    // int       (in): LWORK

    // int (out): INFO
    extern int dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

// Linear algebra functions and utilities
void printMatrix(std::string msg, std::vector<double> vectorV) {
    // Print Vector
    std::cout << msg << " Size=" << vectorV.size() << ", [ ";
    for (double v: vectorV) {
        std::cout << v << " ";
    }
    std::cout << "]" << std::endl;
    return;
}

void printMatrix(std::string msg, std::vector<std::vector<double>> matrixM) {
    // Print Matrix
    std::cout << msg << " Size=" << matrixM.size() << "x" << matrixM[0].size() << ", [ ";
    for (std::vector<double> row: matrixM) {
        std::cout << "[";
        for (double m: row) {
            std::cout << m << " ";
        }
        std::cout << "] ";
    }
    std::cout << "]" << std::endl;
    return;
}

std::vector<double> multMatrix(std::vector<std::vector<double>> matrixM, std::vector<double> vectorV) {
    // Matrix Vector Multiplication
    int sizeMn = matrixM.size(), sizeMm = matrixM[0].size(), sizeVn = vectorV.size();
    // Sanity check
    if (sizeMm != sizeVn) {
        throw std::invalid_argument(std::to_string(sizeMn) + "x" + std::to_string(sizeMm) + " matrix can not be multiplied with a " + std::to_string(sizeVn) + " dimensional vector!");
    }

    std::vector<double> result(sizeMn);
    for (int i=0; i < sizeMn; i++) {
        result[i] = multMatrix(matrixM[i], vectorV);
    }
    return result;
}

double multMatrix(std::vector<double> vectorV, std::vector<double> vectorW) {
    // Vector, Dot Product
    int sizeVn = vectorV.size(), sizeWn = vectorW.size();
    // Sanity check
    if (sizeVn != sizeWn) {
        throw std::invalid_argument("Can not take the dot product of a " + std::to_string(sizeVn) + " dimensional vector with a " + std::to_string(sizeWn) + " dimensional vector!");
    }

    double result = 0;
    for (int i=0; i < sizeVn; i++) {
        result += vectorV[i]*vectorW[i];
    }
    return result;
}

std::vector<double> multMatrix(std::vector<double> vectorV, double scalarS) {
    // Vector, Scalar Multiplication
    int sizeVn = vectorV.size();
    std::vector<double> result(sizeVn);
    for (int i=0; i < sizeVn; i++) {
        result[i] = scalarS*vectorV[i];
    }
    return result;
}

std::vector<double> addMatrix(std::vector<double> vectorV, std::vector<double> vectorW) {
    // Vector Addition
    int sizeVn = vectorV.size(), sizeWn = vectorW.size();
    // Sanity check
    if (sizeVn != sizeWn) {
        throw std::invalid_argument("Can not add a " + std::to_string(sizeVn) + " dimensional vector with a " + std::to_string(sizeWn) + " dimensional vector!");
    }

    std::vector<double> result(sizeVn);
    for (int i=0; i < sizeVn; i++) {
        result[i] = vectorV[i] + vectorW[i];
    }
    return result;
}

// Spectrum Analysis: Eigenvalues, Fourier
void eigenvalues(std::vector<std::vector<double>> matrixA) {
    int sizeMn = matrixA.size(), sizeMm = matrixA[0].size();
    // Sanity check
    if (sizeMn != sizeMm) {
        throw std::invalid_argument("The " + std::to_string(sizeMn) + "x" + std::to_string(sizeMm) + " matrix is not square. Cannot find the eigenvalues!");
    }

    // Convert to a contiguous array
    double data[sizeMn*sizeMn];
    for (int i = 0; i < sizeMn; i++) {
        std::memcpy(data + i*sizeMn, &matrixA[i][0], sizeMn*sizeof(double));
    }

    // Setup parameters for DGEEV
    char Nchar='N';
    std::vector<double> eigReal(sizeMn, 0);
    std::vector<double> eigImag(sizeMn, 0);
    int one=1;
    int lwork=3*sizeMn;
    double* work=new double[lwork];
    int info;

    dgeev_(&Nchar, &Nchar, &sizeMn, data, &sizeMn, eigReal.data(), eigImag.data(), nullptr, &one, nullptr, &one, work, &lwork, &info);

    if (info != 0) {
        std::cout << "DGEEV: Failed with error code -> " << info << std::endl;
    }

    //Check eigenvalues
    //the eigenvalues should be real negative numbers.
    for (double vi: eigImag) {
        if (vi != 0) {
            throw std::invalid_argument("Error: Complex eigenvalue found in system matrix! All eigenvalues should be real negative numbers.");
        }
    }

    for (double vr: eigReal) {
        if (vr > 0) {
            throw std::invalid_argument("Error: Positive eigenvalue found in system matrix! All eigenvalues should be real negative numbers.");
        }
    }

    //Frequency Distribution
    std::vector<double> frequency(sizeMn, 0);
    for (int i = 0; i < sizeMn; i++) {
        frequency[i] = std::sqrt(std::abs(eigReal[i]))/(2*M_PI);
    }
    plotHistogram(frequency, std::max(10, sizeMn/10));
    plotHistogram(frequency, (int)std::sqrt(sizeMn));
}

void powerSpectrum(std::vector<double> t, std::vector<std::vector<double>> data) {
    int sizeN = t.size();
    int numMass = data[0].size();

    std::vector<double> f(sizeN, 0);
    std::vector<double> powerDensity(sizeN, 0);

    double dt = t[1] - t[0];

    // Initialize f
    for (int i = 0; i < sizeN; i++) {
        f[i] = i/(sizeN*dt);
    }

    // Setup FFTW
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sizeN);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sizeN);

    p = fftw_plan_dft_1d(sizeN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int k = 1; k < numMass - 1; k++) {
        double meanValue = 0;
        for (std::vector<double> values: data) {
            meanValue += values[k];
        }
        meanValue /= sizeN;

        for (int i = 0; i < sizeN; i++) {
            in[i][0] = data[i][k] - meanValue; // Remove the DC offset
            in[i][1] = 0;
        }

        fftw_execute(p); // Take the FFT of the contents pointed to by `in`

        for (int i = 0; i < sizeN/2; i++) {
            powerDensity[i] += (out[i][0]*out[i][0] + out[i][1]*out[i][1])/(sizeN*dt)/(numMass - 2);
        }
        std::cout << "\rFourier Transform: " << k << "/" << numMass - 2;
        std::cout.flush();
    }
    std::cout << std::endl;

    // Show Input data
    showPlot(f, powerDensity, "Power Spectrum", "Frequency f", "Power Spectral Density");

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

// Plotting and Graphing
void showPlot(std::vector<double> t, std::vector<double> y, std::string title, std::string xlabel, std::string ylabel, std::string graphlabel, double xMin, double xMax, double yMin, double yMax) {
    Gnuplot gp;
    gp << "set title '" << title << "'" << std::endl;
    gp << "set xlabel '" << xlabel << "'" << std::endl;
    gp << "set ylabel '" << ylabel << "'" << std::endl;

    if (xMin != xMax) {
        gp << "set xrange [" << xMin << ":" << xMax << "]" << std::endl;
    }

    if (yMin != yMax) {
        gp << "set yrange [" << yMin << ":" << yMax << "]" << std::endl;
    }

    if (graphlabel == "") {
        graphlabel = "notitle";
    } else {
        graphlabel = "title '" + graphlabel + "'";
    }

    gp << "plot" << gp.file1d(boost::make_tuple(t, y)) << "with lines " << graphlabel << std::endl;
    gp << "pause mouse close" << std::endl;
}

void plotHistogram(std::vector<double> vectorV, int bins) {
    int sizeVn = vectorV.size();
    std::sort(vectorV.begin(), vectorV.end());

    //Find bin edges
    double minV = vectorV[0], maxV = vectorV.back();
    double dV = (maxV - minV)/(bins);
    std::vector<double> bin_edges(bins + 1, 0);
    for (int i = 0; i < bins + 1; i++) {
        bin_edges[i] = minV + i*dV;
    }

    //Count elements
    std::vector<double> counts(bins, 0);
    double binMin, binMax;

    int root = -1;
    int oldPivot = root;
    int pivot = (sizeVn - root)/2 + root;

    for (int i = 0; i < bins; i++) {
        binMin = bin_edges[i];
        binMax = bin_edges[i + 1];

        int loop_fallback = 0;
        while (loop_fallback < 2*sizeVn) {
            //Check if pivot straddles the interval
            if (((binMax < vectorV[pivot + 1]) || (pivot == sizeVn)) && ((vectorV[pivot] <= binMax) || (pivot == -1))) {
                break;
            } else if (((vectorV[pivot + 1] <= binMax) || (pivot == sizeVn)) && ((vectorV[pivot] <= binMax) || (pivot == -1))) {
                oldPivot = pivot;
                pivot = (sizeVn - pivot)/2 + pivot;
            } else if (((binMax < vectorV[pivot + 1]) || (pivot == sizeVn)) && ((binMax < vectorV[pivot]) || (pivot == -1))) {
                pivot = (pivot - oldPivot)/2 + oldPivot;
            }

            //Check if pivot is maxed out
            if ((pivot == oldPivot) && (pivot == sizeVn - 1)) {
                break;
            }
            loop_fallback++;
        }

        if (loop_fallback > sizeVn) {
            throw std::invalid_argument("Histogram: maximum loop count exceeded!");
        }
        counts[i] = pivot - root;
        root = pivot;
        oldPivot = root;
        pivot = (sizeVn - root)/2 + root;
    }

    //Define histogram curve
    double yMin = 0, yMax = *max_element(counts.begin(), counts.end());
    double xMin = bin_edges.front(), xMax = bin_edges.back();
    double graphMargin = 0.1;

    double xMargin = graphMargin*(xMax - xMin)/2;
    double yMargin = graphMargin*(yMax - yMin)/2;

    std::vector<double> xValues(2*bins + 4, 0);
    std::vector<double> yValues(2*bins + 4, 0);

    xValues[0] = xMin - xMargin;
    xValues[1] = xMin;
    xValues[2*bins + 2] = xMax;
    xValues[2*bins + 3] = xMax + xMargin;
    for (int i = 0; i < bins; i++) {
        xValues[2*i + 2] = bin_edges[i];
        xValues[2*i + 3] = bin_edges[i + 1];
        yValues[2*i + 2] = counts[i];
        yValues[2*i + 3] = counts[i];
    }

    //Graph the histogram using Gnuplot
    showPlot(xValues, yValues, "Frequency Distribution", "Frequency f", "Counts", "", xMin - xMargin, xMax + xMargin, yMin - yMargin, yMax + yMargin);
    return;
}

void renderFrame(std::string fname, double t, std::vector<double> x, double energy) {
    // Graph the springs using Gnuplot
    Gnuplot gp;
    gp << "set terminal png" << std::endl;
    gp << "set output '" + fname + "'" << std::endl;
    gp << boost::format("set title 'Spring Evolution: Time %.2f'") % t << std::endl;
    gp << "set yrange [-1:1]" << std::endl;
    gp << "set xrange [" << x[0] << ":" << x[x.size() - 1] <<"]" << std::endl;
    gp << "set xlabel 'Position x'" << std::endl;
    gp << "plot" << gp.file1d(boost::make_tuple(x, std::vector<double>(x.size(), 0))) << boost::format(" with points pointtype 4 title 'Spring System: Energy %.4e',\\") % energy << std::endl;
    gp <<           gp.file1d(boost::make_tuple(x, std::vector<double>(x.size(), 0))) << " with lines dashtype '--' notitle" << std::endl;
    return;
}

void animation(std::string fname, std::vector<double> t, std::vector<std::vector<double>> x, std::vector<double> energy, int stride) {
    using namespace std::filesystem;

    int sizeT = t.size();
    std::string tmp_dir = "tmp";

    // Create temp directory if it does not exist
    if (!exists(tmp_dir)) {
        create_directory(tmp_dir);
    }

    int width = std::to_string(sizeT - 1).size(); // get fix with for index
    for (int i = 0; stride*i < sizeT; i++) {
        std::string frame_name = std::to_string(stride*i);
        frame_name = fname + "_" + std::string(width - frame_name.size(), '0') + frame_name + ".png";

        renderFrame(tmp_dir + "/" + frame_name, t[stride*i], x[stride*i], energy[stride*i]);
        std::cout << "\rRender Frames: " << stride*(i + 1) << "/" << sizeT;
        std::cout.flush();
    }
    std::cout << std::endl;
}

// Simulation
void simulation(std::string fname, double m, double k, std::vector<double>& x_0, std::vector<double>& x_dot_0, double tMax, int steps, int verbose) {
    // Sanity check
    int num_mass = x_0.size();
    if (x_0.size() != x_dot_0.size()) {
        throw std::invalid_argument("Vectors x_0 and x_dot_0 must have the same length!");
    }

    // Initialize Dynamic Arrays
    double dt = tMax/(steps - 1);

    std::vector<double> t(steps);
    for (int i = 0; i < steps; i++) {
        t[i] = dt*i;
    }

    std::vector<std::vector<double>> x(steps, x_0);
    std::vector<std::vector<double>> x_dot(steps, x_dot_0);
    std::vector<double> energy(steps, 0.0);
    energy[0] = 0.5*m*multMatrix(x_dot[0], x_dot[0]) + 0.5*k*multMatrix(x[0], x[0]);

    // Initialize System Configuration Array
    std::vector<std::vector<double>> matrixA(num_mass, std::vector<double>(num_mass));
    for (int i = 1; i < num_mass - 1; i++) {
        matrixA[i][i - 1] = k/m;
        matrixA[i][i] = -2*k/m;
        matrixA[i][i + 1] = k/m;
    }

    // DEBUG
    eigenvalues(matrixA);

    // Simulation
    for (int i = 1; i < steps; i++) {
        x_dot[i] = addMatrix(x_dot[i - 1], multMatrix(multMatrix(matrixA, x[i - 1]), dt));
        x[i] = addMatrix(x[i - 1], multMatrix(x_dot[i], dt));

        energy[i] = 0.5*m*multMatrix(x_dot[i], x_dot[i]) + 0.5*k*multMatrix(x[i], x[i]);

        std::cout << "\rSimulation step: " << i << "/" << steps - 1;
        std::cout.flush();
    }
    std::cout << std::endl;

    // Graph the power spectrum and the energy using Gnuplot
    powerSpectrum(t, x);
    showPlot(t, energy, "Energy Evolution", "Time t", "Energy", "energy");

    animation(fname, t, x, energy, 10);
    return;
}

int main() {
    int num_mass = 50, steps = 10000;
    double m = 1.0, k = 1.0, tMax = 500.0;

    // Initialize initial values
    std::vector<double> x_0(num_mass);
    for (int i = 0; i < num_mass; i++) {
        x_0[i] = 1.0*i;
    }
    std::vector<double> x_dot_0(num_mass, 0.0);

    x_dot_0[1] = 3;

    simulation("1D_spring_system", m, k, x_0, x_dot_0, tMax, steps);
    return 0;
}
