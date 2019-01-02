#ifndef FUNCTIONS_H
#define FUNCTIONS_H

int gaussian_random(unsigned int n, double *random_num, double sigma=1, double avg=0);
int uniform_random(unsigned int n, double *random_num, double r_min, double r_max);
int gaussian_random_adjust(unsigned int n, double *random_num, double sigma, double avg=0);
int uniform_random_adjust(unsigned int n, double *random_num, double avg=0);
bool iszero(double &x);
double rd(double x, double y, double z);
#endif // FUNCTIONS_H
