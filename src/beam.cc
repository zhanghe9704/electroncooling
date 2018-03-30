
#include "beam.h"
#include <cmath>
#include <cstring>

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle): charge_number_(charge_number), mass_number_(mass_number),
           kinetic_energy_(kinetic_energy), emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), sigma_s_(sigma_s),
           particle_number_(n_particle) {
    mass_ = mass_number*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = (sigma_s_>0)?true:false;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
}

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double n_particle): charge_number_(charge_number), mass_number_(mass_number), kinetic_energy_(kinetic_energy),
           emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), particle_number_(n_particle) {
    mass_ = mass_number_*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = false;
    sigma_s_ = -1;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
}

int Beam::set_center(int i, double x) {
    if(i<3) {
        center_[i] = x;
        return 0;
    }
    else {
        perror("Error index for electron beam center!");
        return 1;
    }
}

int UniformCylinder::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    for(unsigned int i=0; i<n_particle; ++i){
        if(x[i]*x[i]+y[i]*y[i]<=r2) ne[i] = density;
    }
    return 0;
}

int UniformCylinder::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                             double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    for(unsigned int i=0; i<n_particle; ++i){
        if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2) ne[i] = density;
    }
    return 0;
}

int UniformHollow::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle) {
    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));

    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}

int UniformHollow::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                             double cy, double cz) {
    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));
     //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}



int UniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>= left_end && x[i]*x[i]+y[i]*y[i]<=r2)
            ne[i] = density;
    }
    return 0;
}

int UniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                          double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*ebeam.beta()*k_c);

    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end && (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2)
            ne[i] = density;
    }
    return 0;
}

int UniformHollowBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    memset(ne, 0, n_particle*sizeof(double));

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(z[i]<=right_end && z[i]>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}

int UniformHollowBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                          double cy, double cz){

    int nq = ebeam.charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*ebeam.beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        double z_shifted = z[i]+cz;
        if(z_shifted<=right_end && z_shifted>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
    return 0;
}



int EllipticUniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle){

    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*ebeam.beta()*k_c);
    memset(ne, 0, n_particle*sizeof(double));
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>=left_end && inv_rh2*x[i]*x[i]+inv_rv2*y[i]*y[i]<=1)
            ne[i] = density;
    }
    return 0;
}

int EllipticUniformBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle,
                                  double cx, double cy, double cz){
    int nq = ebeam.charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*ebeam.beta()*k_c);

    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    memset(ne, 0, n_particle*sizeof(double));
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(unsigned int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end &&
           inv_rh2*(x[i]+cx)*(x[i]+cx)+inv_rv2*(y[i]+cy)*(y[i]+cy)<=1)
            ne[i] = density;
    }
    return 0;
}

EBeam::EBeam(double gamma, double tmp_tr, double tmp_long, EBeamShape &shape_defined):
    Beam(-1, k_me/k_u, (gamma-1)*k_me, 0, 0, 0, 0),tmp_tr_(tmp_tr),tmp_long_(tmp_long){
    shape_ = &shape_defined;
    v_rms_long_ = sqrt(tmp_long_/this->mass())*0.001*k_c;
    v_rms_tr_ = sqrt(tmp_tr/this->mass())*0.001*k_c;
}

int GaussianBunch::density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n_particle){
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_y_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    for(unsigned int i=0; i<n_particle; ++i){
//        ne[i] = amp*exp(-0.5*(x[i]*x[i]/sigma_x2+y[i]*y[i]/sigma_y2+z[i]*z[i]/sigma_s2));
        ne[i] = amp*exp(x[i]*x[i]*sigma_x2+y[i]*y[i]*sigma_y2+z[i]*z[i]*sigma_s2);
    }
    return 0;
}

int GaussianBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx,
                           double cy, double cz){
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_y_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    //ion_center - electron_center
    cx -= ebeam.center(0);
    cy -= ebeam.center(1);
    cz -= ebeam.center(2);
    for(unsigned int i=0; i<n_particle; ++i){
        ne[i] = amp*exp((x[i]+cx)*(x[i]+cx)*sigma_x2+(y[i]+cy)*(y[i]+cy)*sigma_y2+(z[i]+cz)*(z[i]+cz)*sigma_s2);
//        ne[i] = amp*exp(-0.5*((x[i]+cx)*(x[i]+cx)/sigma_x2+(y[i]+cy)*(y[i]+cy)/sigma_y2+(z[i]+cz)*(z[i]+cz)/sigma_s2));
    }
    return 0;
}
