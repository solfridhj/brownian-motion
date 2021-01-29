#include <iostream>
#include <armadillo>
#include "Particle.h"

using namespace arma;
using namespace std;


int main(int argc, char *argv[]) {
    // Seeding random functions. Use the same value for reproducibility

    double t_end;
    double x_start;
    double deltaU;
    double delta_t;

    t_end = 1;

    x_start = 0.0;

    delta_t = 0.2*pow(10, -4);
    //---------------------------------------------------------------------------

    // Convert to SI-units
    x_start = x_start*pow(10, -6); // [m]

    // Set other parameters
    double tau = 0.78;


    double alpha = 0.2;

    double L = 20*pow(10, -6); // [m] length

    double k_BT = 26*pow(10, -3)*1.602176565*pow(10, -19);  // [J]

    //deltaU = 0.1*k_BT; // [J]

    deltaU = 80*1.602176565*pow(10, -19); // [J];


    double eta = pow(10, -3); // [Pa s]
    double ri = 3*12*pow(10, -9); // [m]

    // sets flashing to be true
    bool flashing = true;

    Particle p(alpha, L, tau, eta, ri, deltaU, k_BT, delta_t, x_start, t_end, flashing);

    arma::arma_rng::set_seed(40);

    // Runs simulation, saves the output as CSV
    //p.eulerScheme();
    //p.plotUr();
    //p.drift();
    //p.particleDistribution();
    p.ensemble();

}
