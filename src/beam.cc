
#include "beam.h"
#include <cmath>
#include <cstring>
#include "arbitrary_electron_beam.h"

#include <fstream>

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
    energy_spread_ = beta_*beta_*dp_p_;
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
    energy_spread_ = beta_*beta_*dp_p_;
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
    v_rms_tr_ = sqrt(tmp_tr_/this->mass())*0.001*k_c;
    if(shape_defined.shape()==Shape::PARTICLE_BUNCH) {
        velocity_ = Velocity::USER_DEFINE;
        temperature_ = Temperature::USER_DEFINE;
    }
}

EBeam::EBeam(double gamma, EBeamShape &shape_defined):Beam(-1, k_me/k_u, (gamma-1)*k_me, 0, 0, 0, 0) {
    shape_ = &shape_defined;
    tmp_tr_ = 0;
    tmp_long_ = 0;
    v_rms_long_ = sqrt(tmp_long_/this->mass())*0.001*k_c;
    v_rms_tr_ = sqrt(tmp_tr_/this->mass())*0.001*k_c;
    if(shape_defined.shape()==Shape::PARTICLE_BUNCH) {
        velocity_ = Velocity::USER_DEFINE;
        temperature_ = Temperature::USER_DEFINE;
    }
}

int GaussianBunch::density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n_particle){
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_x_);
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
    double sigma_x2 = -1/(2*sigma_x_*sigma_x_);
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

//ParticleBunch::ParticleBunch(double n_electron, std::string filename, unsigned long int n, double length, int line_skip,
//                             bool binary, int buffer, int s, double neutralisation):n_electron_(n_electron),
//                             filename_(filename), n_(n), length_(length), line_skip_(line_skip), binary_(binary),
//                             buffer_(buffer), s_(s), neutralisation_(neutralisation){
////    x.reserve(n);
////    y.reserve(n);
////    z.reserve(n);
////    vx.reserve(n);
////    vy.reserve(n);
////    vz.reserve(n);
//
//    n_ = load_electrons(x, y, z, vx, vy, vz, filename, n, line_skip, binary, buffer);
////    list_e_.reserve(n_);
//    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
//}

//ParticleBunch::ParticleBunch(double n_electron, std::string filename, unsigned long int n, int line_skip, bool binary,
//                             int buffer, int s, double neutralisation): n_electron_(n_electron), filename_(filename), n_(n),
//                             line_skip_(line_skip), binary_(binary), buffer_(buffer), s_(s), neutralisation_(neutralisation){
////    x.reserve(n);
////    y.reserve(n);
////    z.reserve(n);
////    vx.reserve(n);
////    vy.reserve(n);
////    vz.reserve(n);
//
//    n_ = load_electrons(x, y, z, vx, vy, vz, filename, n, line_skip, binary, buffer);
////    list_e_.reserve(n_);
//
//    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
//    auto itr = z.begin();
//    double z_max = *itr;
//    double z_min = *itr;
//    ++itr;
//    for(; itr!=z.end(); ++itr) {
//        if(*itr>z_max) z_max = *itr;
//        if(*itr<z_min) z_min = *itr;
//    }
//    length_ = z_max - z_min;
//}
//
//
//ParticleBunch::ParticleBunch(double n_electron, std::string filename, double length, int line_skip,
//                             bool binary, int buffer, int s, double neutralisation):n_electron_(n_electron),
//                             filename_(filename), n_(0), length_(length), line_skip_(line_skip), binary_(binary),
//                             buffer_(buffer), s_(s), neutralisation_(neutralisation){
////    x.reserve(n);
////    y.reserve(n);
////    z.reserve(n);
////    vx.reserve(n);
////    vy.reserve(n);
////    vz.reserve(n);
//
//    n_ = load_electrons(x, y, z, vx, vy, vz, filename, 0, line_skip, binary, buffer);
////    list_e_.reserve(n_);
//    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
//}
//
//ParticleBunch::ParticleBunch(double n_electron, std::string filename, int line_skip, bool binary,
//                             int buffer, int s, double neutralisation): n_electron_(n_electron), filename_(filename), n_(0),
//                             line_skip_(line_skip), binary_(binary), buffer_(buffer), s_(s), neutralisation_(neutralisation){
////    x.reserve(n);
////    y.reserve(n);
////    z.reserve(n);
////    vx.reserve(n);
////    vy.reserve(n);
////    vz.reserve(n);
//
//    n_ = load_electrons(x, y, z, vx, vy, vz, filename, 0, line_skip, binary, buffer);
////    list_e_.reserve(n_);
//
//    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
//    auto itr = z.begin();
//    double z_max = *itr;
//    double z_min = *itr;
//    ++itr;
//    for(; itr!=z.end(); ++itr) {
//        if(*itr>z_max) z_max = *itr;
//        if(*itr<z_min) z_min = *itr;
//    }
//    length_ = z_max - z_min;
//}

void ParticleBunch::load_particle(long int n) {
    n_ = load_electrons(x, y, z, vx, vy, vz, filename_, n, line_skip_, binary_, buffer_);
    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
    if(length_==0) {
        auto itr = z.begin();
        double z_max = *itr;
        double z_min = *itr;
        ++itr;
        for(; itr!=z.end(); ++itr) {
            if(*itr>z_max) z_max = *itr;
            if(*itr<z_min) z_min = *itr;
        }
        length_ = z_max - z_min;
    }
}

void ParticleBunch::load_particle() {
    load_particle(0);
}

int ParticleBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n) {
    double rate = n_electron_/n_;
    std::vector<unsigned int> list_i;
    unsigned int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_)
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_z, v_rms_t, v_rms_l);}
    else
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);}
//    tpr_l.clear();
//    tpr_t.clear();
//    double c2_inv = 1.0/(k_c*k_c);
//    for(auto& v: v_rms_l) tpr_l.push_back(k_me*v*v*1e6*c2_inv);
//    for(auto& v: v_rms_t) tpr_t.push_back(k_me*v*v*1e6*c2_inv);

    for(unsigned int i=0; i<n; ++i) ne[i] *= rate;

    return 0;
}

int ParticleBunch::density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy,
                       double cz) {
    double rate = n_electron_/n_;
    for(unsigned int i=0; i<n; ++i) {
        x[i] -= cx;
        y[i] -= cy;
        z[i] -= cz;
    }
    std::vector<unsigned int> list_i;
    unsigned int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_)
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_z, v_rms_t, v_rms_l);}
    else
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);}
    for(unsigned int i=0; i<n; ++i) {
        x[i] += cx;
        y[i] += cy;
        z[i] += cz;
    }
    double c2_inv = 1.0/(k_c*k_c);
    for(auto& v: v_rms_l) tpr_l.push_back(k_me*v*v*1e6*c2_inv);
    for(auto& v: v_rms_t) tpr_t.push_back(k_me*v*v*1e6*c2_inv);
    for(unsigned int i=0; i<n; ++i) ne[i] *= rate;
    return 0;
}
