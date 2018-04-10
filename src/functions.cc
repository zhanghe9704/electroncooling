#include <chrono>
#include <random>

int gaussian_random(unsigned int n, double *random_num, double sigma, double avg){
    std::default_random_engine generator;
//    generator.seed(0);
    generator.seed(rand());
    std::normal_distribution<double> distribution(avg,sigma);
    for(unsigned int i=0; i<n; ++i) random_num[i] = distribution(generator);
    return 0;
}

int uniform_random(unsigned int n, double *random_num, double r_min, double r_max){
    std::default_random_engine generator;
//    generator.seed(0);
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


