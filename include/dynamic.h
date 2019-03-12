#ifndef DYNAMIC_H
#define DYNAMIC_H

#include <fstream>
#include <iostream>
#include <string>
#include "beam.h"
#include "ecooling.h"
#include "ibs.h"
#include "ring.h"
#include "cooler.h"

using std::string;
enum class DynamicModel {RMS, PARTICLE, MODEL_BEAM = PARTICLE, TURN_BY_TURN};

class Luminosity {
    double dx_ = 0;
    double dy_ = 0;
    double np_1_ = 0;
    double np_2_ = 0;
    double freq_ = 1;
    double sigma_x1_ = 0;
    double sigma_y1_ = 0;
    double sigma_x2_ = 0;
    double sigma_y2_ = 0;
    double bet_x1_ = 0;
    double bet_y1_ = 0;
    double bet_x2_ = 0;
    double bet_y2_ = 0;
    double geo_emit_x1_ = 0;
    double geo_emit_y1_ = 0;
    double geo_emit_x2_ = 0;
    double geo_emit_y2_ = 0;
    bool use_ion_emittance_ = true;
public:
    void set_distance(double dx, double dy){dx_=dx; dy_=dy;}
    void set_freq(double f){freq_=f;}
    bool set_use_ion_emit(bool b){use_ion_emittance_ = b;}
    void set_geo_emit(double emit_x, double emit_y, int i);
    void set_beam_size(double sigma_x, double sigma_y, int i);
    void set_particle_number(double n, int i);
    void set_bet(double bet_x, double bet_y, int i);
    bool use_ion_emittance(){return use_ion_emittance_;}
    double luminosity();
};

class DynamicParas{
    double time_;
    int n_step_;
    int n_sample_;
    double dt_;
    bool ibs_ = true;
    bool ecool_ = true;
    bool fixed_bunch_length_ = false;
    bool reset_time_ = true;
    bool overwrite_ = true;
    bool calc_luminosity_ = false;
    int output_intvl_ = 1;
    int ion_save_intvl_ = -1;
    string filename_ = "output_dynamic.txt";
    DynamicModel model_ = DynamicModel::RMS;
//    int n_sample_;
 public:
    Twiss twiss_ref;
    int n_step(){return n_step_;}
    double time(){return time_;}
    double dt(){return dt_;}
    bool ibs(){return ibs_;}
    bool ecool(){return ecool_;}
    bool fixed_bunch_length(){return fixed_bunch_length_;}
    bool reset_time(){return reset_time_;}
    bool overwrite(){return overwrite_;}
    bool calc_lum(){return calc_luminosity_;}
    int output_intval(){return output_intvl_;}
    int ion_save_intvl(){return ion_save_intvl_;}
    int n_sample(){return n_sample_;}
//    int n_sample() {assert(model_==DynamicModel::MODEL_BEAM); return n_sample_;};
    DynamicModel model(){return model_;}
    int set_model(DynamicModel model){model_ = model; return 0;}
    int set_ion_save(int x){ion_save_intvl_ = x; return 0;}
    int set_output_file(string filename){filename_ = filename; return 0;}
    int set_output_intvl(int x){output_intvl_ = x; return 0;}
    int set_n_sample(int x){n_sample_ = x; return 0;}
    int set_fixed_bunch_length(bool b){fixed_bunch_length_ = b; return 0;}
    int set_reset_time(bool b){reset_time_ = b; return 0;}
    int set_overwrite(bool b) {overwrite_ = b; return 0;}
    int set_calc_lum(bool b) {calc_luminosity_ = b;}
    string output_file(){return filename_;}
    DynamicParas(double time, int n_step):time_(time),n_step_(n_step){dt_ = time_/n_step_;}
    DynamicParas(double time, int n_step, bool ibs, bool ecool):
        time_(time),n_step_(n_step),ibs_(ibs),ecool_(ecool){dt_ = time_/n_step_;}
    DynamicParas(double time, int n_step, bool ibs, bool ecool, DynamicModel model):
        time_(time),n_step_(n_step),ibs_(ibs),ecool_(ecool), model_(model){dt_ = time_/n_step_;}
};

//int dynamic(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring, std::ofstream &outfile);
int dynamic(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring);
#endif // DYNAMIC_H
