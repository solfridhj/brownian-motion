//
// Created by Solfrid Johansen on 27/02/2020.
//

#ifndef ASSIGNMENT_2_PARTICLE_H
#define ASSIGNMENT_2_PARTICLE_H

#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;


class Particle {

private:
    double alpha;   // Asymmetry factor, between 0 and 1
    double L;       // [m], Period of the sawtooth potential
    double tau;     // [s], Period of flashing
    double gamma;   // Gamma, drag coefficient
    double deltaU;  // [V] Strength of factor
    double k_BT;       // [K] Temperature in kelvin
    double delta_t;
    double delta_t_hat;
    double eta;
    double r_i;

    double xHat_n;
    double tHat_n;

    double omega;
    double DHat;

    // End position
    double t_end;
    double t_hat_end;

    // Start position
    double x_start;

    // Number of simulations to average ofer
    int nr_trails;

    // data
    vec t;
    vec x;
    vec x_average;
    vec x_end;

    bool flashing;

    double f(double t);     // Flashing potential (time dep. part)
    double Ur(double x);    // Asymmetric saw-tooth (position dep. part)
    double Fr(double x);    // Position-dependent part of the force

    bool checkCriterion();


public:

    Particle(double alpha, double L, double tau, double nu, double r_i, double deltaU, double k_BT, double delta_t, double x_start, double t_end, bool flashing);
    double U(double x, double t);       // Returns the reduced potential at x and t

    void plotUr();      // Saves as csv
    void plotf();     // Saves as csv
    void plotFr();      // Saves as csv

    double F(double x, double t);       // Value of force, acts in the x direction

    static double randomNumber();       // Generates a random number

    // For producing many random numbers to ensure correct behaviour
    double generateRandomNumbers();

    double euler(double n, double x_n);

    void eulerScheme();

    // Computes particle distribution over several realizations
    void particleDistribution();


    // Computes average drift velocity
    void computeDriftVel();

    void drift();

    void ensemble();

};


#endif //ASSIGNMENT_2_PARTICLE_H
