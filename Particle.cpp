//
// Created by Solfrid Johansen on 27/02/2020.
//

#include "Particle.h"



Particle::Particle(double alpha, double L, double tau, double eta, double r_i, double deltaU, double k_BT, double delta_t, double x_start, double t_end, bool flashing): x_start{x_start}, t_end{t_end}, alpha{alpha}, L{L}, delta_t{delta_t}, flashing{flashing}, eta{eta}, r_i{r_i}{
    // Converting all to SI units

    // Value of gamma
    gamma = 6*datum::pi*eta*r_i;


    this->t_end = t_end;
    DHat = k_BT/(deltaU);
    this->k_BT = k_BT;
    this->deltaU = deltaU;
    // Computing omega
    omega = deltaU/(gamma*L*L);

    // mulig noe reduced med dette...
    this-> tau = tau;
    // Reduced end-time
    this->t_end = 100;
    t_hat_end = t_end*omega;
    t_hat_end = 100;
    // Reduced time-step
    delta_t_hat = delta_t*omega;

}


double Particle::Ur(double x){
    // Compute U_r in reduced units
    x = abs(fmod(x, 1));

    if( (x >= 0 &&  x < alpha)){
        return x/alpha;
    }else if((x >= alpha && x < 1)){
        return (1 - x) / (1 - alpha);
    }else{
        cerr << "Error in computing the Ur potential. " << endl;
    }

    return -100;
}

void Particle::plotUr(){

    vec x_vec = linspace(-3, 3, 200);
    vec U_vec(200);

    for(int i = 0; i < x_vec.n_rows; i++){
        U_vec(i) = Ur(x_vec(i));
    }

    U_vec.save("../data/potential/potential_data.csv", csv_ascii);
    x_vec.save("../data/potential/x_data.csv", csv_ascii);
}


double Particle::f(double t){

    t = fmod(t, tau*omega);

    if( (t >= 0) && (t < 3.0*tau/4.0*omega)){

        return 0;

    }else if((t >= 3.0*tau/4.0*omega) && (t<=tau*omega)){

        return 1;

    }
    cerr << "Error in computing f" << endl;
    return -100;

}


void Particle::plotf(){
    vec x = linspace(0, 40, 200);
    vec flash(200);

    for(int i = 0; i < x.n_rows; i++){
        flash(i) = f(x(i));
    }

    ofstream pot{"../data/potential/flashing_data.csv"};
    for (int i = 0; i < x.n_rows; i++) {
        pot << flash(i) << ",";
    }

    ofstream lin{"../data/potential/x_flashing_data.csv"};
    for (int i = 0; i < x.n_rows; i++) {
        lin << x(i) << ",";
    }
}

double Particle::U(double x, double t){
    // Returns the potential energy as a function  of x and t
    return Ur(x)*f(t);
}


double Particle::Fr(double x){

    // Modulates the value
    x = abs(fmod(x, 1));

    // Value of force depends on position of the particle
    if( (x >= 0 &&  x < alpha)){
        return -1.0/alpha;
    }else if((x >= alpha && x<1)){
        return (1.0)/(1.0 - alpha);
    }else{

        cerr << "Error in computing the force for " << x << endl;
    }
    return 0;
}


void Particle::plotFr(){

    vec x = linspace(0, 5, 200);
    vec force(200);

    for(int i = 0; i < x.n_rows; i++){
        force(i) = Fr(x(i));
    }

    ofstream pot{"../data/potential/force_data.csv"};
    for (int i = 0; i < x.n_rows; i++) {
        pot << force(i) << ",";
    }

    ofstream lin{"../data/potential/x_force_data.csv"};
    for (int i = 0; i < x.n_rows; i++) {
        lin << x(i) << ",";
    }
}

double Particle::F(double x, double t){
    // Include flashing

    if(flashing){
        return Fr(x)*f(t);
    }else{
        return Fr(x);
    }
}

double Particle::randomNumber(){
    return randn<double>();
}

double Particle::generateRandomNumbers(){
    int n = 1000000; // Nr samples
    // Store samples
    vec samples(n);

    for(int i = 0; i < n; i++){
        samples(i) = randomNumber();
    }

    samples.save("../data/random/rand_vals.csv", csv_ascii);
}

bool Particle::checkCriterion(){
    // Computes a good value of delta_t

    double dUdx{0};
    double crit{0};

    // Finds the maximum value of the force, for a given alpha
    if(1.0/alpha > 1.0/(1.0-alpha)){
        dUdx = 1.0/alpha;
    }else{
        dUdx = 1.0/(1.0-alpha);
    }

    crit = dUdx*delta_t_hat + 4*sqrt(2*DHat*delta_t_hat);
    if (crit < 0.1*alpha) {

        return true;
    } else {
        cout << crit << " but should be < 0.1*alpha: " << 0.1*alpha << endl;
        //return false;
        return true;
    }

}

double Particle::euler(double n, double x_n){
    // Compute the value of the next x, n+1
    double epsilon = randomNumber();
    // SET TO 1 FOR POTENTIAL, 0 FOR NO POTENTIAL
    int potential = 1;
    // Reduced time
    tHat_n = n*delta_t*omega;
    double x_n_1 = 0;

    if(potential == 1){
        x_n_1 = x_n + F(x_n, tHat_n)*delta_t_hat + sqrt(2*DHat*delta_t_hat)*epsilon;
    }else{
        // NB; this is NOT in reduced units, due to deltaU = 0 not allowing this.
        x_n_1 = x_n + sqrt(2*k_BT*delta_t/gamma)*epsilon;
    }
    return x_n_1;
}

void Particle::eulerScheme(){
    // Make it so what the actual values are converted to reduced units,
    // then computed and returned as actual values

    // In seconds
    double t_start = 0;

    if(!checkCriterion()){
        cout << "Too high value of delta_t, select smaller " << endl;
        return;
    }
    // Compute average over many eulerSchemes
    cout << "Computing Euler Scheme" << endl;
    // Array of time values to evaluate the position at. Needed for data

    // For storing the average over all computations
    x_average = regspace(t_start, delta_t_hat, t_end);


    //t = regspace(t_start, delta_t_hat, t_end);
    t = regspace(t_start, delta_t_hat, t_hat_end);
    // First value is in reduced units

    x = regspace(t_start, delta_t_hat, t_end);
    x(0) = x_start/L;
    for(int i = 0; i < x.n_rows-1; i++){
        x(i+1) = euler(i, x(i));
        //cout << i  << " " << x.n_rows -1<< endl;
    }

    // UNCOMMENT TO COMPUTE AVERAGES
    /*
    nr_trails = 10;
    // To store the end position of each trail
    x_end.set_size(nr_trails);
    for(int k = 0; k < nr_trails; k++){

        // Store last element
        x_end(k) = x(x.n_rows-1);
        x_average += x;
        cout << k << endl;
    }

    x_average = x_average / nr_trails*L;

    x_average.save("../data/single_particle/x_mean.csv", csv_ascii);
     */

    // Save data as CSV

    int tau_ = static_cast<int>(tau);
    ofstream xVal{"../data/single_particle/x_val.csv"};
    for (int i = 0; i < x.n_rows; i++) {
        //xVal << x(i)*L<< ",";
        xVal << x(i)<< ",";
    }

    ofstream tVal{"../data/single_particle/t_val.csv"};
    for (int i = 0; i < x.n_rows; i++) {

        tVal << t(i) << ",";
    }

}

void Particle::particleDistribution(){
    // Compute the average of the distribution of energy

    eulerScheme();
    vec U(nr_trails);
    double scaling  = deltaU/L;
    // should reach some equilibrium first
    for(int i = 0; i < nr_trails; i++){
        // U is a value between 0 and 1.
        U(i) = 1-abs(Ur(x_end(i)));
    }

    U.save("../data/boltzmann/U_avg_10.csv", csv_ascii);
}


void Particle::drift() {

    double tau_min = 0.4;
    double tau_max = 6;
    double steps = 0.1;

    vec taus = regspace(tau_min, steps, tau_max);
    vec vel(taus.n_rows);

    int count = 0;
    for (double i = tau_min; i < tau_max; i = i + steps) {
        tau = i;
        // 1000 realizations
        eulerScheme();
        // Find velocity by average end position/time
        vel(count) = x_average(x_average.n_rows - 1)/t_hat_end;
        cout << count << endl;
        count ++;
    }
    taus.save("../data/single_particle/taus_r2.csv", csv_ascii);
    vel.save("../data/single_particle/vel_r2.csv", csv_ascii);
}

void Particle::ensemble(){

    // Ensemble size
    int N = 300;

    eulerScheme();
    // Matrix to save all computed states
    mat x_ensemble(x.n_rows, N);

    for(int i = 0; i < N; i++){
        // Computes eulerScheme for one particle
        eulerScheme();
        // Store the computed position value as a row (in metres)
        x_ensemble.col(i) = x;
    }
    x_ensemble.save("../data/many_particles/ensemble_r2.csv", csv_ascii);
    t.save("../data/many_particles/time_r2.csv", csv_ascii);
}