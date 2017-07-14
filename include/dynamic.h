#ifndef DYNAMIC_H
#define DYNAMIC_H

#include <fstream>
#include <iostream>
#include "beam.h"
#include "ecooling.h"
#include "ibs.h"
#include "ring.h"
#include "cooler.h"

enum class DynamicModel {RMS, MODEL_BEAM, TURN_BY_TURN};

class DynamicParas{
    double time_;
    int n_step_;
    double dt_;
    bool ibs_ = true;
    bool ecool_ = true;
    DynamicModel model_ = DynamicModel::RMS;
//    int n_sample_;
 public:
    int n_step(){return n_step_;}
    double time(){return time_;}
    double dt(){return dt_;}
    bool ibs(){return ibs_;}
    bool ecool(){return ecool_;}
//    int n_sample() {assert(model_==DynamicModel::MODEL_BEAM); return n_sample_;};
    DynamicModel model(){return model_;}
    int set_model(DynamicModel model){model_ = model; return 0;}
    DynamicParas(double time, int n_step):time_(time),n_step_(n_step){dt_ = time_/n_step_;}
    DynamicParas(double time, int n_step, bool ibs, bool ecool):
        time_(time),n_step_(n_step),ibs_(ibs),ecool_(ecool){dt_ = time_/n_step_;}
    DynamicParas(double time, int n_step, bool ibs, bool ecool, DynamicModel model):
        time_(time),n_step_(n_step),ibs_(ibs),ecool_(ecool), model_(model){dt_ = time_/n_step_;}
//    DynamicParas(double time, int n_step, bool ibs, bool ecool, int n_sample):
//        time_(time),n_step_(n_step),ibs_(ibs),ecool_(ecool), model_(DynamicModel::MODEL_BEAM), n_sample_(n_sample)
//        {dt_ = time_/n_step_;}
};

int dynamic(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring, std::ofstream &outfile);

#endif // DYNAMIC_H
