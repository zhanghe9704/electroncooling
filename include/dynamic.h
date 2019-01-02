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

class DynamicParas{
    double time_;
    int n_step_;
    int n_sample_;
    double dt_;
    bool ibs_ = true;
    bool ecool_ = true;
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
