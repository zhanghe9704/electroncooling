#include "ring.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

Lattice::Lattice(std::string filename) {
    std::ifstream infile;
    infile.open(filename.c_str());
    if(!infile) std::cout<<"Error: failed to load the lattice file!"<<std::endl;
    std::string line;

    bool read = false;
    unsigned long int cnt = 0;
    while(std::getline(infile,line)){
        std::istringstream iss(line);
        std::string name;
        std::string keyword;
        double num;
        iss>>name;
//        if (name.compare("\"MACHINE$START\"")==0) read=true;
        std::size_t found = name.find("$START");
        if (found!=std::string::npos) read=true;

        if(read){
            iss>>keyword;
            iss>>num;
            if(s_.empty()||num>s_.at(cnt-1)){
                s_.push_back(num);
                iss>>num;
                betx_.push_back(num);
                iss>>num;
                alfx_.push_back(num);
                iss>>num;
                mux_.push_back(num);
                iss>>num;
                dx_.push_back(num);
                iss>>num;
                dpx_.push_back(num);
                iss>>num;
                bety_.push_back(num);
                iss>>num;
                alfy_.push_back(num);
                iss>>num;
                muy_.push_back(num);
                iss>>num;
                dy_.push_back(num);
                iss>>num;
                dpy_.push_back(num);
                ++cnt;
            }
        }
        found = name.find("$END");
        if (found!=std::string::npos){
//        if(name.compare("\"MACHINE$END\"")==0){
            n_element_ = cnt;
            read = false;
            break;
        }
    }
    circ_ = s_.at(n_element_-1) - s_.at(0);
    for(int i=0; i<n_element_-1; ++i){
        l_element_.push_back(s_.at(i+1)-s_.at(i));
    }
    l_element_.push_back(l_element_.at(0));
}

Ring::Ring(double circ, Beam &beam_defined):circ_(circ) {
    beam_ = &beam_defined;
    lattice_ = nullptr;
    if(beam_->bunched())
        beta_s_ = beam_->sigma_s()/beam_->dp_p();
}

Ring::Ring(Lattice &lattice_defined, Beam &beam_defined) {
    beam_ = &beam_defined;
    lattice_ = &lattice_defined;
    circ_ = lattice_->circ();
    if(beam_->bunched())
        beta_s_ = beam_->sigma_s()/beam_->dp_p();
}


//
//Beam::Beam(int n, double A, double KE, double emit_nx, double emit_ny, double dp_p, double sigma_s, double n_particle):n_charge(n),A(A),KE(KE),emit_nx(emit_nx),emit_ny(emit_ny),dp_p(dp_p),sigma_s(sigma_s),n_particle(n_particle){
//    mass = A*u;
//    gamma = 1+KE/mass;
//    beta = sqrt(gamma*gamma-1)/gamma;
//    r = ke*n_charge*n_charge*e*1e-6/mass;
//    bunched = (sigma_s>0)?true:false;
//    emit_x = emit_nx/(beta*gamma);
//    emit_y = emit_ny/(beta*gamma);
//
//}
//
//Beam::Beam(int n, double A, double KE, double emit_nx, double emit_ny, double dp_p, double n_particle):n_charge(n),A(A),KE(KE),emit_nx(emit_nx),emit_ny(emit_ny),dp_p(dp_p),n_particle(n_particle){
//    mass = A*u;
//    gamma = 1+KE/mass;
//    beta = sqrt(gamma*gamma-1)/gamma;
//    r = ke*n_charge*n_charge*e*1e-6/mass;
//    bunched = false;
//    sigma_s = -1;
//    emit_x = emit_nx/(beta*gamma);
//    emit_y = emit_ny/(beta*gamma);
//
//}
//
//
////Ring::Ring(double gamma_t, int h, double V, std::string filename, Beam &beam_defined):gamma_t(gamma_t),h(h),V(V){
////    lattice.lattice_from_file(filename);
////    beam = &beam_defined;
////
////    double gamma = beam->get_gamma();
////    eta = 1/(gamma_t*gamma_t) - 1/(gamma*gamma);
////
////    double Q = h*e*fabs(eta)*beam->get_Z()*V/(2*pi*beam->get_mass_J()*beam->get_A()*gamma);
////    Qs = sqrt(Q)/beam->get_beta();
////    Bs = lattice.get_circ()*fabs(eta)/(2*pi*Qs);
////
////}
//
//Ring::Ring(double gamma_t, int h, double V, double circ, Beam &beam_defined):gamma_t(gamma_t),h(h),V(V),circ(circ){
//    beam = &beam_defined;
//
//    double gamma = beam->get_gamma();
//    eta = 1/(gamma_t*gamma_t) - 1/(gamma*gamma);
//
//    double Q = h*e*fabs(eta)*beam->get_Z()*V/(2*pi*beam->get_mass_J()*beam->get_A()*gamma);
//    Qs = sqrt(Q)/beam->get_beta();
//    Bs = circ*fabs(eta)/(2*pi*Qs);
//    beam->set_sigma_s(Bs*beam->get_dp_p());
//
//}
//
//Ring::Ring(double gamma_t, int h, double V, Lattice &lattice_defined, Beam &beam_defined):gamma_t(gamma_t),h(h),V(V){
//
//    lattice = &lattice_defined;
//    beam = &beam_defined;
//    circ = lattice->get_circ();
//
//    double gamma = beam->get_gamma();
//    eta = 1/(gamma_t*gamma_t) - 1/(gamma*gamma);
//
//    double Q = h*e*fabs(eta)*beam->get_Z()*V/(2*pi*beam->get_mass_J()*beam->get_A()*gamma);
//    Qs = sqrt(Q)/beam->get_beta();
//    Bs = circ*fabs(eta)/(2*pi*Qs);
//    beam->set_sigma_s(Bs*beam->get_dp_p());
//
//}
//
//Cooler::Cooler(double length, double n_section, double B, double beta_h, double beta_v, double disp_h, double disp_v, double alpha_h, double alpha_v, double der_disp_h, double der_disp_v)
//    :length(length),n_section(n_section),B(B),beta_h(beta_h),beta_v(beta_v),disp_h(disp_h),disp_v(disp_v),alpha_h(alpha_h),alpha_v(alpha_v),der_disp_h(der_disp_h),der_disp_v(der_disp_v){
//    }
//
//
//double UniformCylinder::get_density(double x, double y, double z, Beam &beam){
//
//    double density = 0;
//    if(x*x+y*y<=radius*radius){
//        int nq = beam.get_n_charge();
//        if (nq<0) nq *= -1;
//
//        density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//    }
//    return density;
//}
//
//
//int UniformCylinder::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if(x[i]*x[i]+y[i]*y[i]<=radius*radius) ne[i] = density;
//    }
//    return 0;
//}
//
//int UniformCylinder::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n, double cx, double cy, double cz){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    //ion_center - electron_center
//    cx -= beam.get_center(0);
//    cy -= beam.get_center(1);
//    cz -= beam.get_center(2);
//
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius) ne[i] = density;
//    }
//    return 0;
//}
//
//
//
//double UniformCylinderBunch::get_density(double x, double y, double z, Beam &beam){
//    double density = 0;
//    if (z<0.5*length&&z>-0.5*length) {
//        if(x*x+y*y<=radius*radius) {
//            int nq = beam.get_n_charge();
//            if (nq<0) nq *= -1;
//
//            density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//        }
//    }
//    return density;
//}
//
//
//int UniformCylinderBunch::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if(z[i]<0.5*length&&z[i]>-0.5*length)
//            if(x[i]*x[i]+y[i]*y[i]<=radius*radius) ne[i] = density;
//    }
//    return 0;
//}
//
//int UniformCylinderBunch::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n, double cx, double cy, double cz){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    //ion_center - electron_center
//    cx -= beam.get_center(0);
//    cy -= beam.get_center(1);
//    cz -= beam.get_center(2);
//
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if((z[i]+cz)<0.5*length&&(z[i]+cz)>-0.5*length)
//            if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius) ne[i] = density;
//    }
//    return 0;
//}
//
//
//
//double darkCurrent::get_density(double x, double y, double z, Beam &beam){
//    double density = 0;
//    if (z<0.5*peak_length&&z>-0.5*peak_length) {
//        if(x*x+y*y<=radius*radius) {
//            int nq = beam.get_n_charge();
//            if (nq<0) nq *= -1;
//            density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//        }
//    }
//    else if(z<0.5*(length) || z>-0.5*(length)) {
//        if(x*x+y*y<=radius*radius) {
//            int nq = beam.get_n_charge();
//            if (nq<0) nq *= -1;
//            density = dark_current/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//        }
//    }
//    return density;
//}
//
//
//int darkCurrent::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//    double density_dark = dark_current/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if(z[i]<0.5*peak_length&&z[i]>-0.5*peak_length) {
//            if(x[i]*x[i]+y[i]*y[i]<=radius*radius) ne[i] = density;
//        }
//        else if (z[i]<0.5*(length) || z[i]>-0.5*(length)) {
//             if(x[i]*x[i]+y[i]*y[i]<=radius*radius) ne[i] = density_dark;
//        }
//    }
//    return 0;
//}
//
//int darkCurrent::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n, double cx, double cy, double cz){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//    double density_dark = dark_current/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    //ion_center - electron_center
//    cx -= beam.get_center(0);
//    cy -= beam.get_center(1);
//    cz -= beam.get_center(2);
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if((z[i]+cz)<0.5*peak_length&&(z[i]+cz)>-0.5*peak_length) {
//            if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius) ne[i] = density;
//        }
//        else if (z[i]+cz<0.5*(length) || z[i]+cz>-0.5*(length)) {
//            if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius) ne[i] = density_dark;
//        }
//    }
//    return 0;
//}
//
//
//
//
//
//double UniformCylinderSlope::get_density(double x, double y, double z, Beam &beam){
//
//    double density = 0;
//    if(x*x+y*y<=radius*radius){
//        int nq = beam.get_n_charge();
//        if (nq<0) nq *= -1;
//
//        density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//    }
//
//    if((z>0.5*length-slope_begin*length)&&(z<=0.5*length+slope_begin*length)){
//        density *= (0.5*length+slope_begin*length-z)/(2*slope_begin*length);
//    }
//    else if((z>=-0.5*length-slope_end*length)&&(z<-0.5*length+slope_end*length)){
//        density *= (z+0.5*length+slope_end*length)/(2*slope_end*length);
//    }
//
//    return density;
//}
//
//int UniformCylinderSlope::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if(x[i]*x[i]+y[i]*y[i]<=radius*radius) ne[i] = density;
//    }
//
//    double length_begin_1 = 0.5*length-slope_begin*length;
//    double length_begin_2 = 0.5*length+slope_begin*length;
//    double slope_length_1 = 2*slope_begin*length;
//    double length_end_1 = -0.5*length-slope_end*length;
//    double length_end_2 = -0.5*length+slope_end*length;
//    double slope_length_2 = 2*slope_end*length;
//    for(unsigned int i=0; i<n; ++i){
//        if((z[i]>length_begin_1)&&(z[i]<=length_begin_2)){
//        ne[i] *= (length_begin_2-z[i])/slope_length_1;
//        }
//        else if((z[i]>=length_end_1)&&(z[i]<length_begin_2)){
//            ne[i] *= (z[i]-length_end_1)/slope_length_2;
//        }
//    }
//    return 0;
//}
//
//int UniformCylinderSlope::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n, double cx, double cy, double cz){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    //ion_center - electron_center
//    cx -= beam.get_center(0);
//    cy -= beam.get_center(1);
//    cz -= beam.get_center(2);
//
//
//    memset(ne, 0, n*sizeof(double));
//    for(unsigned int i=0; i<n; ++i){
//        if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius) ne[i] = density;
//    }
//
//    double length_begin_1 = 0.5*length-slope_begin*length;
//    double length_begin_2 = 0.5*length+slope_begin*length;
//    double slope_length_1 = 2*slope_begin*length;
//    double length_end_1 = -0.5*length-slope_end*length;
//    double length_end_2 = -0.5*length+slope_end*length;
//    double slope_length_2 = 2*slope_end*length;
//    for(unsigned int i=0; i<n; ++i){
//        if(((z[i]+cz)>length_begin_1)&&((z[i]+cz)<=length_begin_2)){
//        ne[i] *= (length_begin_2-z[i]-cz)/slope_length_1;
//        }
//        else if(((z[i]+cz)>=length_end_1)&&((z[i]+cz)<length_begin_2)){
//            ne[i] *= (z[i]+cz-length_end_1)/slope_length_2;
//        }
//    }
//    return 0;
//}
//
//
//double GaussianCylinder::get_density(double x, double y, double z, Beam &beam){
//
//    double density = 0;
//    if((x*x+y*y<=radius*radius)&&(z*z<=9*sigma_s*sigma_s)){
//        int nq = beam.get_n_charge();
//        if (nq<0) nq *= -1;
//
//        density = I_peak/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//        density *= exp(-0.5*z*z/(sigma_s*sigma_s));
//    }
//    return density;
//}
//
//
//int GaussianCylinder::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I_peak/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    memset(ne, 0, n*sizeof(double));
//    double coef_exp = -1/(2*sigma_s*sigma_s);
//    for(unsigned int i=0; i<n; ++i){
//        if((x[i]*x[i]+y[i]*y[i]<=radius*radius)&&(z[i]*z[i]<=9*sigma_s*sigma_s)) ne[i] = density*exp(z[i]*z[i]*coef_exp);
//    }
//    return 0;
//}
//
//int GaussianCylinder::get_density(double *x, double *y, double *z, Beam &beam, double *ne, unsigned int n, double cx, double cy, double cz){
//
//    int nq = beam.get_n_charge();
//    if (nq<0) nq *= -1;
//    double density = I_peak/(pi*radius*radius*nq*e*beam.get_beta()*c0);
//
//    //ion_center - electron_center
//    cx -= beam.get_center(0);
//    cy -= beam.get_center(1);
//    cz -= beam.get_center(2);
//
//
//    memset(ne, 0, n*sizeof(double));
//    double coef_exp = -1/(2*sigma_s*sigma_s);
//    for(unsigned int i=0; i<n; ++i){
//        if(((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=radius*radius)&&((z[i]+cz)*(z[i]+cz)<=9*sigma_s*sigma_s))
//            ne[i] = density*exp(coef_exp*(z[i]+cz)*(z[i]+cz));
//    }
//    return 0;
//}
//
//
////EBeamCooling::EBeamCooling(double gamma, double T_tr, double T_long, Shape &shape_defined):gamma(gamma),T_tr(T_tr),T_long(T_long){
////    shape = &shape_defined;
////    beta = sqrt(1-1/(gamma*gamma));
////    v_rms_long = sqrt(T_long/mass)*0.001*c0;
////    v_rms_tr = sqrt(T_tr/mass)*0.001*c0;
////}
//
//EBeam::EBeam(double gamma, double T_tr, double T_long, Shape &shape_defined):Beam(-1, me/u, (gamma-1)*me, 0, 0, 0, 0),T_tr(T_tr),T_long(T_long){
//    shape = &shape_defined;
//    v_rms_long = sqrt(T_long/this->get_mass())*0.001*c0;
//    v_rms_tr = sqrt(T_tr/this->get_mass())*0.001*c0;
//}
//
//bool EBeam::get_bunched(){
//
//    bool bunched;
//
//    switch(shape->get_shape()){
//    case ShapeList::GaussianBunch:{
//        bunched = true;
//        break;
//    }
//    default:{
//        bunched = false;
//        break;
//    }
//    }
//
//    return bunched;
//
//}
//
//double GaussianBunch::get_density(double x, double y, double z, Beam &ebeam){
//    return n/(sqrt(8*pi*pi*pi)*sigma_x*sigma_y*sigma_s)*exp(-1/2*(x*x/(sigma_x*sigma_x)+y*y/(sigma_y*sigma_y)+z*z/(sigma_s*sigma_s)));
//}
//
//int GaussianBunch::get_density(double *x, double *y, double *z,  Beam &ebeam, double *ne, unsigned int n_particle){
//
//    double A = n/(sqrt(8*pi*pi*pi)*sigma_x*sigma_y*sigma_s);
//    double sigma_x2 = sigma_x*sigma_x;
//    double sigma_y2 = sigma_y*sigma_y;
//    double sigma_s2 = sigma_s*sigma_s;
//    for(unsigned int i=0; i<n_particle; ++i){
//        ne[i] = A*exp(-0.5*(x[i]*x[i]/sigma_x2+y[i]*y[i]/sigma_y2+z[i]*z[i]/sigma_s2));
//    }
//
//    return 0;
//}
//
//int GaussianBunch::get_density(double *x, double *y, double *z,  Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy, double cz){
//
//    double A = n/(sqrt(8*pi*pi*pi)*sigma_x*sigma_y*sigma_s);
//    double sigma_x2 = sigma_x*sigma_x;
//    double sigma_y2 = sigma_y*sigma_y;
//    double sigma_s2 = sigma_s*sigma_s;
//
//    //ion_center - electron_center
//    cx -= ebeam.get_center(0);
//    cy -= ebeam.get_center(1);
//    cz -= ebeam.get_center(2);
//
//    for(unsigned int i=0; i<n_particle; ++i){
//        ne[i] = A*exp(-1/2*((x[i]+cx)*(x[i]+cx)/sigma_x2+(y[i]+cy)*(y[i]+cy)/sigma_y2+(z[i]+cz)*(z[i]+cz)/sigma_s2));
//    }
//
//    return 0;
//}
//

