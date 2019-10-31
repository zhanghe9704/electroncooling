#include <algorithm>
#include <assert.h>
#include <chrono>
#include <math.h>
#include <random>

int gaussian_random(unsigned int n, double *random_num, double sigma, double avg){
    std::default_random_engine generator;
//    generator.seed(2);
    generator.seed(rand());
    std::normal_distribution<double> distribution(avg,sigma);
    for(unsigned int i=0; i<n; ++i) random_num[i] = distribution(generator);
    return 0;
}

int uniform_random(unsigned int n, double *random_num, double r_min, double r_max){
    std::default_random_engine generator;
//    generator.seed(2);
    generator.seed(rand());
    std::uniform_real_distribution<double> uniform_dis(r_min,r_max);
    for(unsigned int i=0; i<n; ++i) random_num[i] = uniform_dis(generator);
    return 0;
}

int gaussian_random_adjust(unsigned int n, double *random_num, double sigma, double avg) {
    double mean = 0;
    for(unsigned int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;

    double sigma_calc = 0;
    for(unsigned int i=0; i<n; ++i) {
        random_num[i] -= mean;
        sigma_calc += random_num[i]*random_num[i];
    }
    sigma_calc = sqrt(sigma_calc/n);

    double adjust_width = sigma/sigma_calc;
    for(unsigned int i=0; i<n; ++i) {
        random_num[i] *= adjust_width;
        random_num[i] += avg;
    }
    return 0;
}

int uniform_random_adjust(unsigned int n, double *random_num, double avg) {
    double mean = 0;
    for(unsigned int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;
    double adjust = avg - mean;
    for(unsigned int i=0; i<n; ++i) random_num[i] += adjust;
    return 0;
}

bool iszero(double &x) {
    double err = 1.0e-60;
    if (x<err&&x>-err) {
//        x = 0;
        return true;
    }
    return false;
}


/** \brief Computes Carlson's elliptic integral of the second kind, RD(x, y, z).
 * x and y must be nonnegative, and at most one can be zero. z must be positive.
 * TINY must be at least twice the negative 2/3 power of the machine overflow limit.
 * BIG must be at most 0.1 ¡Á ERRTOL times the negative 2/3 power of the machine underflow limit.
 * The algorithm follows the example in "Numerical recipes in C: the art of scientific computing" by Press W.H.,
 * Teukolsky S.A., Vetterling W.T., and Flannery B.P.
 *
 * \param x nonnegative, only one of x and y can be zero.
 * \param y nonnegative, only one of x and y can be zero.
 * \param z positive.
 * \return the value of the integral.
 *
 */

double rd(double x, double y, double z) {
    double ERRTOL = 0.01;
    double TINY = 1.0e-25;
    double BIG = 4.5e21;
    double C1 = 3.0/14.0;
    double C2 = 1.0/6.0;
    double C3 = 9.0/22.0;
    double C4 = 3.0/26.0;
    double C5 = 0.25*C3;
    double C6 = 1.5*C4;

    if(std::min(x,y)<0.0 || std::min(x+y,z)<TINY || std::max(std::max(x,y),z)>BIG)
        assert(false&&"Invalid arguments in rd!");
    double sum = 0;
    double fac = 1.0;
    double delta_x, delta_y, delta_z, ave;
    do {
        double sqrt_x = sqrt(x);
        double sqrt_y = sqrt(y);
        double sqrt_z = sqrt(z);
        double alamb = sqrt_x*(sqrt_y+sqrt_z)+sqrt_y*sqrt_z;
        sum += fac/(sqrt_z*(z+alamb));
        fac *= 0.25;
        x=0.25*(x+alamb);
        y=0.25*(y+alamb);
        z=0.25*(z+alamb);
        ave=0.2*(x+y+3.0*z);
        delta_x = (ave-x)/ave;
        delta_y = (ave-y)/ave;
        delta_z = (ave-z)/ave;
    }
    while (std::max(std::max(fabs(delta_x),fabs(delta_y)), fabs(delta_z))>ERRTOL);
    double ea = delta_x*delta_y;
    double eb = delta_z*delta_z;
    double ec = ea - eb;
    double ed = ea - 6.0*eb;
    double ee = ed + ec + ec;
    return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delta_z*ee)+delta_z*(C2*ee+delta_z*(-C3*ec+delta_z*C4*ea)))/(ave*sqrt(ave));
}
