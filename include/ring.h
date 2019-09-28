#ifndef RING_H
#define RING_H

#include <cassert>
#include <memory>
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
    double s(int i) const {return s_.at(i);}
    double betx(int i) const {return betx_.at(i);}
    double alfx(int i) const {return alfx_.at(i);}
    double mux(int i) const {return mux_.at(i);}
    double dx(int i) const {return dx_.at(i);}
    double dpx(int i) const {return dpx_.at(i);}
    double bety(int i) const {return bety_.at(i);}
    double alfy(int i) const {return alfy_.at(i);}
    double muy(int i) const {return muy_.at(i);}
    double dy(int i) const {return dy_.at(i);}
    double dpy(int i) const {return dpy_.at(i);}
    int n_element() const {return n_element_;}
    double l_element(int i) const {return l_element_.at(i);}
    double circ() const {return circ_;}
    Lattice(std::string filename);

};

struct Tunes {
    double qx = 0;
    double qy = 0;
    double qs = 0;
};

struct RF {
    double v = 0; //Voltage in [V].
    int h = 1; //Harmonic number
    double phi = 0; //phase in [2*PI]
    double gamma_tr = 0; //Transition gamma
};

class Ring{
    double beta_s_ = 0;         //Synchrotron function, use to calculate rms bunch length from momentum spread
    double circ_ = 0;        //Circumference of the ring;
    double f0_ = 0;               // revolution frequency.
    double w0_ = 0;       // angular revolution frequency.
    double slip_factor_ = 0;     //slip factor.
 public:
    Beam *beam_;
    Lattice *lattice_;
    Tunes tunes;
    RF rf;
//    Tunes *tunes = nullptr;
//    RF *rf = nullptr;
//    std::shared_ptr<Beam> beam_;
//    std::shared_ptr<Lattice> lattice_;
    double beta_s(){assert(beam_->bunched()); return beta_s_;}
    double circ(){return circ_;}
    double f0(){return f0_;}
    double w0(){return w0_;}
    double slip_factor(){return slip_factor_;}
    double calc_rf_voltage();
    double calc_sync_tune_by_rf();
    Ring(double circ, Beam &beam_defined);
    Ring(Lattice &lattice_defined, Beam &beam_defined);
    void set_rf();
};

struct Twiss{
    double bet_x = 0;
    double bet_y = 0;
    double alf_x = 0;
    double alf_y = 0;
    double disp_x = 0;
    double disp_y = 0;
    double disp_dx = 0;
    double disp_dy = 0;
};

#endif // RING_H
