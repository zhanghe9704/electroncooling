#ifndef BEAM_H
#define BEAM_H

#include "constants.h"
#include <cstdio>
#include <memory>

class Beam{
    int charge_number_;   //Number of charges
    double mass_number_;       //mass = A*u [MeV/c^2]
    double mass_;    //unit in MeV/c^2
    double r_;       //classical radius, in m
    double kinetic_energy_;      //kinetic energy, in MeV
    double beta_;    //Lorentz factors
    double gamma_;   //Lorentz factors
    double emit_nx_; //normalized horizontal emittance, in m
    double emit_ny_; //normalized vertical emittance, in m
    double emit_x_;  //geometrical horizontal emittance, in m
    double emit_y_;  //geometrical vertical emittance, in m
    double dp_p_;     //momentum spread dp/p
    double sigma_s_; //RMS bunch length. set it to -1 for coasting beam, in m
    double particle_number_; //number of particles
    bool bunched_;   //Return true if beam is bunched.
    double center_[3] = {0,0,0};

public:
    int set_emit_nx(double x){emit_nx_ = x; emit_x_ = emit_nx_/(beta_*gamma_); return 0;}
    int set_emit_ny(double x){emit_ny_ = x; emit_y_ = emit_ny_/(beta_*gamma_); return 0;}
    int set_emit_x(double x){emit_x_ = x; emit_nx_ = beta_*gamma_*emit_x_; return 0;}
    int set_emit_y(double x){emit_y_ = x; emit_ny_ = beta_*gamma_*emit_y_; return 0;}
    int set_dp_p(double x){dp_p_ = x; return 0;}
    int set_sigma_s(double x){sigma_s_ = x; return 0;}
    int set_center(double cx, double cy, double cz){center_[0] = cx; center_[1] = cy; center_[2] = cz; return 0;}
    int set_center(int i, double x);
    int charge_number(){return charge_number_;}
    double mass(){return mass_;}
    double kinetic_energy(){return kinetic_energy_;}
    double beta(){return beta_;}
    double gamma(){return gamma_;}
    double emit_nx(){return emit_nx_;}
    double emit_ny(){return emit_ny_;}
    double emit_x(){return emit_x_;}
    double emit_y(){return emit_y_;}
    double dp_p(){return dp_p_;}
    double sigma_s(){return sigma_s_;}
    double r(){return r_;}
    double particle_number(){return particle_number_;}
    double mass_number(){return mass_number_;}
    double mass_J(){return mass_*1e6*k_e;}
    bool bunched(){return bunched_;}
    int center(double &cx, double &cy, double &cz){cx = center_[0]; cy = center_[1]; cz = center_[2]; return 0;}
    double center(int i){ if (i<3) return center_[i]; else perror("Error index for electron beam center!"); return 1.0;}
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double sigma_s, double n_particle);
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double n_particle);
};

enum class Shape {UNIFORM_CYLINDER, GAUSSIAN_BUNCH, UNIFORM_BUNCH, GAUSSIAN_CYLINDER, ELLIPTIC_UNIFORM_BUNCH,
    UNIFORM_HOLLOW, UNIFORM_HOLLOW_BUNCH};

enum class Velocity {CONST, USER_DEFINE, SPACE_CHARGE}  ;

class EBeamShape{
 public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    virtual int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n)=0;
    //(cx, cy, cz) is the relative position of the ion beam center to the electron beam center.
    virtual int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx,
                            double cy, double cz)=0;
    virtual Shape shape()=0;
    virtual bool bunched() = 0;
    virtual double length()=0; //For bunched electron beam, return full length of the electron bunch.
    virtual double neutralisation()=0;
    EBeamShape(){};
};

class UniformCylinder: public EBeamShape{
    double current_;                   //Current of the beam in A
    double radius_;              //Radius of the beam in meter
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double radius(){return radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_CYLINDER;}
    double length(){perror("length() not defined for UniformCylinder, which is coasting"); return 0;}
    bool bunched(){return false;}
    UniformCylinder(double current, double radius, double neutralisation=2):current_(current),radius_(radius),
                    neutralisation_(neutralisation){};
};

class UniformHollow: public EBeamShape {
    double current_;    //Peak current, the current as if the beam is coasting.
    double in_radius_;
    double out_radius_;
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double out_radius(){return out_radius_;}
    double in_radius(){return in_radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_HOLLOW;}
    double length(){perror("length() not defined for UniformHollow, which is coasting"); return 0;}
    bool bunched(){return false;}
    UniformHollow(double current, double in_radius, double out_radius, double neutralisation=2):current_(current),
        in_radius_(in_radius), out_radius_(out_radius),neutralisation_(neutralisation){};
};

class UniformHollowBunch: public EBeamShape {
    double current_;
    double in_radius_;
    double out_radius_;
    double neutralisation_;
    double length_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    double current(){return current_;}
    double out_radius(){return out_radius_;}
    double in_radius(){return in_radius_;}
    double neutralisation(){return neutralisation_;}
    Shape shape(){return Shape::UNIFORM_HOLLOW_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    UniformHollowBunch(double current, double in_radius, double out_radius, double length, double neutralisation=2):current_(current),
        in_radius_(in_radius), out_radius_(out_radius),length_(length), neutralisation_(neutralisation){};
};

class GaussianBunch: public EBeamShape{
    double n_electron_;
    double sigma_x_;
    double sigma_y_;
    double sigma_s_;
    double neutralisation_;
 public:
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n_particle, double cx, double cy,
                double cz);
    Shape shape(){return Shape::GAUSSIAN_BUNCH;}
    double length(){return 6*sigma_s_;}
    bool bunched(){return true;}
    double neutralisation(){return neutralisation_;}
    GaussianBunch(double n_electron, double sigma_x, double sigma_y, double sigma_s):n_electron_(n_electron),
                sigma_x_(sigma_x),sigma_y_(sigma_y),sigma_s_(sigma_s){};

};


class UniformBunch: public EBeamShape{
    double current_;                   //Current of the beam in A, assuming the beam is DC.
    double radius_;              //Radius of the beam in meter
    double length_;
    double neutralisation_;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
    Shape shape(){return Shape::UNIFORM_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    double current(){return current_;}
    double radius(){return radius_;}
    double neutralisation(){return neutralisation_;}
//    UniformCylinder(double I, double radius, double neutralisation):Shape(ShapeList::uniformCylinder),I(I),radius(radius),neutralisation(neutralisation){};
    UniformBunch(double current, double radius, double length, double neutralisation=2):current_(current),radius_(radius),
            length_(length), neutralisation_(neutralisation){};

};

class EllipticUniformBunch: public EBeamShape{
    double current_;
    double rh_;         //half horizontal axis
    double rv_;         //half vertical axis
    double length_;     //bunch length
    double neutralisation_;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n);
    int density(double *x, double *y, double *z, Beam &ebeam, double *ne, unsigned int n, double cx, double cy, double cz);
    Shape shape(){return Shape::ELLIPTIC_UNIFORM_BUNCH;}
    double length(){return length_;}
    bool bunched(){return true;}
    double neutralisation(){return neutralisation_;}
    EllipticUniformBunch(double current, double rh, double rv, double length, double neutralisation=2):current_(current),
            rh_(rh),rv_(rv),length_(length),neutralisation_(neutralisation){};
};

class EBeam:public Beam{
    double tmp_tr_;            //Transverse temperature, in eV
    double tmp_long_;          //Longitudinal temperature, in eV
    double v_rms_tr_;        //Transverse RMS velocity, in m/s
    double v_rms_long_;      //Longitudinal RMS velocity, in m/s
    Velocity velocity_ = Velocity::CONST;

 public:
    EBeamShape *shape_;            //Shape of the electron beam
    double tmp_tr(){return tmp_tr_;}
    double tmp_long(){return tmp_long_;}
    double v_rms_tr(){return v_rms_tr_;}
    double v_rms_long(){return v_rms_long_;}
    int set_velocity(Velocity velocity){velocity_ = velocity; return 0;}
    Velocity velocity(){return velocity_;}

    int emit_nx(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_ny(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_x(){perror("This function is not defined for cooling electron beam"); return 1;}
    int emit_y(){perror("This function is not defined for cooling electron beam"); return 1;}
    int dp_p(){perror("This function is not defined for cooling electron beam"); return 1;}
    int sigma_s(){perror("This function is not defined for cooling electron beam"); return 1;}
    int n_particle(){perror("This function is not defined for cooling electron beam"); return 1;}
    bool bunched(){return shape_->bunched();}
    double length(){return shape_->length();}

    EBeam(double gamma, double tmp_tr, double tmp_long, EBeamShape &shape_defined);
};
#endif // BEAM_H
