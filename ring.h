#ifndef RING_H
#define RING_H

#include <cassert>
#include <string>
#include <vector>
#include "beam.h"
#include "constants.h"

using std::vector;

//enum class ShapeList {uniformCylinder=1, GaussianBunch=2, GaussianCylinder=3, uniformCylinderSlope=4,
//                      uniformCylinderBunch, darkCurrent};

class Lattice{
    vector<double> s_;
    vector<double> betx_;
    vector<double> alfx_;
    vector<double> mux_;
    vector<double> dx_;
    vector<double> dpx_;
    vector<double> bety_;
    vector<double> alfy_;
    vector<double> muy_;
    vector<double> dy_;
    vector<double> dpy_;
    vector<double> l_element_;
    int n_element_;
    double circ_;
 public:
    double s(int i){return s_.at(i);}
    double betx(int i){return betx_.at(i);}
    double alfx(int i){return alfx_.at(i);}
    double mux(int i){return mux_.at(i);}
    double dx(int i){return dx_.at(i);}
    double dpx(int i){return dpx_.at(i);}
    double bety(int i){return bety_.at(i);}
    double alfy(int i){return alfy_.at(i);}
    double muy(int i){return muy_.at(i);}
    double dy(int i){return dy_.at(i);}
    double dpy(int i){return dpy_.at(i);}
    double n_element(){return n_element_;}
    double l_element(int i){return l_element_.at(i);}
    double circ(){return circ_;}
    Lattice(std::string filename);

};

class Ring{
    double beta_s_ = 0;         //Synchrotron function, use to calculate rms bunch length from momentum spread
    double circ_ = 0;        //Circumference of the ring;
 public:
    Beam *beam_;
    Lattice *lattice_;
    double beta_s(){assert(beam_->bunched()); return beta_s_;}
    double circ(){return circ_;}
    Ring(double circ, Beam &beam_defined);
    Ring(Lattice &lattice_defined, Beam &beam_defined);
};


//
//class Beam{
//    int n_charge;   //Number of charges
//    double A;       //mass = A*u [MeV/c^2]
//    double mass;    //unit in MeV/c^2
//    double KE;      //kinetic energy, in MeV
//    double beta;    //Lorentz factors
//    double gamma;   //Lorentz factors
//    double emit_nx; //normalized horizontal emittance, in m
//    double emit_ny; //normalized vertical emittance, in m
//    double emit_x;  //geometrical horizontal emittance, in m
//    double emit_y;  //geometrical vertical emittance, in m
//    double dp_p;     //momentum spread dp/p
//    double sigma_s; //RMS bunch length. set it to -1 for coasting beam, in m
//    double n_particle; //number of particles
//    double r;       //classical radius, in m
//    bool bunched;   //Return true if beam is bunched.
//    double center[3] = {0,0,0};
//
////    int set_n_charge(int n){n_charge = n; return 0;}
////    int set_mass(double m){mass = m; return 0;}
////    int set_gamma(double g){gamma = g; return 0;}
////    int set_beta(double b){beta = b; return 0;}
////    int set_r(double radius){r = radius; return 0;}
////    int set_bunched(bool b){bunched = b; return 0;}
//
//public:
//    int set_emit_nx(double x){emit_nx = x; emit_x = emit_nx/(beta*gamma); return 0;}
//    int set_emit_ny(double x){emit_ny = x; emit_y = emit_ny/(beta*gamma); return 0;}
//    int set_emit_x(double x){emit_x = x; emit_nx = beta*gamma*emit_x; return 0;}
//    int set_emit_y(double x){emit_y = x; emit_ny = beta*gamma*emit_y; return 0;}
//    int set_dp_p(double x){dp_p = x; return 0;}
//    int set_sigma_s(double x){sigma_s = x; return 0;}
//    int get_n_charge(){return n_charge;}
//    double get_mass(){return mass;}
//    double get_KE(){return KE;}
//    double get_beta(){return beta;}
//    double get_gamma(){return gamma;}
//    double get_emit_nx(){return emit_nx;}
//    double get_emit_ny(){return emit_ny;}
//    double get_emit_x(){return emit_x;}
//    double get_emit_y(){return emit_y;}
//    double get_dp_p(){return dp_p;}
//    double get_sigma_s(){return sigma_s;}
//    double get_r(){return r;}
//    double get_n_particle(){return n_particle;}
//    double get_A(){return A;}
//    int get_Z(){return n_charge;}
//    double get_mass_J(){return mass*1e6*e;}
//    bool get_bunched(){return bunched;}
//    int set_center(double cx, double cy, double cz){center[0] = cx; center[1] = cy; center[2] = cz; return 0;}
//    int get_center(double &cx, double &cy, double &cz){cx = center[0]; cy = center[1]; cz = center[2]; return 0;}
//    double get_center(int i){ if (i<3) return center[i]; else perror("Error index for electron beam center!"); return 1.0;}
//    Beam(int n, double A, double KE, double emit_nx, double emit_ny, double dp_p, double sigma_s, double n_particles);
//    Beam(int n, double A, double KE, double emit_nx, double emit_ny, double dp_p, double n_particles);
//
//
//};
//
//class Ring{
//    double eta;         // slip factor or off momentum factor
//    double gamma_t;     //transition energy
//    int h;              //harmonic number
//    double V;           //RF voltage, in V
//    double Qs;          //Synchrotron tune
//    double Bs;         //Synchrotron function, use to calculate rms bunch length from momentum spread
//    double circ;        //Circumference of the ring;
//
//public:
//    Beam *beam;
//    Lattice *lattice;
//    double get_eta(){return eta;}
//    double get_gamma_t(){return gamma_t;}
//    int get_h(){return h;}
//    double get_V(){return V;}
//    double get_Qs(){return Qs;}
//    double get_Bs(){return Bs;}
//    double get_circ(){return circ;}
////    Ring(double gamma_t, int h, double V, std::string filename, Beam &beam_defined);
//    Ring(double gamma_t, int h, double V, double circ, Beam &beam_defined);
//    Ring(double gamma_t, int h, double V, Lattice &lattice_defined, Beam &beam_defined);
//};
//
//class Shape{
////     ShapeList shape;
//
// public:
//    virtual double get_density(double x, double y, double z, Beam &ebeam)=0;
//    virtual int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n)=0;
//    virtual int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz)=0;
//    virtual ShapeList get_shape()=0;
////   virtual bool get_bunched()=0;
////    virtual ShapeList get_shape(){return shape;}
////    Shape(ShapeList shape):shape(shape){};
//    Shape(){};
//};
//
//class UniformCylinder: public Shape{
//    double I;                   //Current of the beam in A
//    double radius;              //Radius of the beam in meter
//    double neutralisation;
//public:
//    //Calculate the charge density for a given position (x,y,z) in Lab frame.
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::uniformCylinder;}
////    bool get_bunched(){return false;}
////    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
//    UniformCylinder(double I, double radius, double neutralisation):I(I),radius(radius),neutralisation(neutralisation){};
//
//};
//
//class UniformCylinderBunch: public Shape{
//    double I;                   //Current of the beam in A
//    double radius;              //Radius of the beam in meter
//    double length;
//    double neutralisation;
//public:
//    //Calculate the charge density for a given position (x,y,z) in Lab frame.
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::uniformCylinderBunch;}
////    bool get_bunched(){return false;}
////    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
//    UniformCylinderBunch(double I, double radius, double length, double neutralisation):I(I),radius(radius),length(length), neutralisation(neutralisation){};
//
//};
//
//class darkCurrent: public Shape{
//    double I;                   //Current of the beam in A
//    double radius;    //Radius of the beam in meter
//    double peak_length;
//    double length;
//    double dark_current;
//    double dark_length;
//    double neutralisation;
//public:
//    //Calculate the charge density for a given position (x,y,z) in Lab frame.
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::uniformCylinderBunch;}
////    bool get_bunched(){return false;}
////    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
//    darkCurrent(double I, double radius, double peak_length, double dark_current, double dark_length, double neutralisation):
//        I(I),radius(radius),peak_length(peak_length), dark_current(dark_current), dark_length(dark_length), neutralisation(neutralisation){
//        length = peak_length+dark_length;}
//
//};
//
//
//class UniformCylinderSlope: public Shape{
//    double I;                   //Current of the beam in A
//    double radius;              //Radius of the beam in meter
//    double neutralisation;
//    double length;
//    double slope_begin;               //percentage of the length
//    double slope_end;
//public:
//    //Calculate the charge density for a given position (x,y,z) in Lab frame.
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::uniformCylinderSlope;}
////    bool get_bunched(){return false;}
////    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
//    UniformCylinderSlope(double I, double radius, double neutralisation, double length, double slope_begin, double slope_end)
//        :I(I),radius(radius),neutralisation(neutralisation),length(length),slope_begin(slope_begin),slope_end(slope_end){};
//
//};
//
//class GaussianBunch: public Shape{
//    double n;             //Number of particles
//    double sigma_x;
//    double sigma_y;
//    double sigma_s;
//public:
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::GaussianBunch;}
////    bool get_bunched(){return true;}
//    GaussianBunch(double n, double sigma_x, double sigma_y, double sigma_s):n(n),sigma_x(sigma_x),sigma_y(sigma_y),sigma_s(sigma_s){};
//};
//
//class GaussianCylinder: public Shape{
//    double I_peak;
//    double radius;
//    double sigma_s;
//public:
//    double get_density(double x, double y, double z, Beam &ebeam);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
//    int get_density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
//    ShapeList get_shape(){return ShapeList::GaussianCylinder;}
//    GaussianCylinder(double I_peak, double radius, double sigma_s):I_peak(I_peak),radius(radius),sigma_s(sigma_s){};
//
//};
//
//class EBeam:public Beam{
//
//    double T_tr;            //Transverse temperature, in eV
//    double T_long;          //Longitudinal temperature, in eV
//    double v_rms_tr;        //Transverse RMS velocity, in m/s
//    double v_rms_long;      //Longitudinal RMS velocity, in m/s
//
//public:
//    Shape *shape;            //Shape of the electron beam
//    double get_T_tr(){return T_tr;}
//    double get_T_long(){return T_long;}
//    double get_v_rms_tr(){return v_rms_tr;}
//    double get_v_rms_long(){return v_rms_long;}
//
//    int get_emit_nx(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_emit_ny(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_emit_x(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_emit_y(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_dp_p(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_sigma_s(){perror("This function is not defined for cooling electron beam"); return 1;}
//    int get_n_particle(){perror("This function is not defined for cooling electron beam"); return 1;}
//    bool get_bunched();
//
//    EBeam(double gamma, double T_tr, double T_long, Shape &shape_defined);
//
//
//};
//
//
//class Cooler{
//    double length;      // in meter
//    double n_section;
//    double B;           // in Tesla
//    double beta_h;      // in meter
//    double beta_v;      // in meter
//    double alpha_h;
//    double alpha_v;
//    double disp_h;
//    double disp_v;
//    double der_disp_h;
//    double der_disp_v;
//
//public:
//    double get_length(){return length;}
//    double get_n_section(){return n_section;}
//    double get_B(){return B;}
//    double get_beta_h(){return beta_h;}
//    double get_beta_v(){return beta_v;}
//    double get_alpha_h(){return alpha_h;}
//    double get_alpha_v(){return alpha_v;}
//    double get_disp_h(){return disp_h;}
//    double get_disp_v(){return disp_v;}
//    double get_der_disp_h(){return der_disp_h;}
//    double get_der_disp_v(){return der_disp_v;}
//    int set_disp_h(double h){disp_h = h; return 0;}
//    int set_disp_v(double v){disp_v = v; return 0;}
//    Cooler(double length, double n_section, double B, double beta_h, double beta_v, double disp_h=0, double disp_v=0, double alpha_h=0, double alpha_v=0, double der_disp_h=0, double der_disp_v=0);
//
//};
#endif // RING_H
