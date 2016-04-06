#ifndef ECOOLING_H
#define ECOOLING_H

#include "beam.h"
#include "cooler.h"
#include "force.h"
#include "ring.h"

enum class IonSample {SINGLE_PARTICLE, MONTE_CARLO};
//enum class ECoolScenario {DC_COAST, DC_BUNCH, BUNCH_COAST, BUNCH_BUNCH};

class EcoolRateParas{
    IonSample ion_sample_ = IonSample::SINGLE_PARTICLE;
    unsigned int n_sample_ = 125000;
    unsigned int n_tr_ = 50;
    unsigned int n_long_ = 50;
    bool shift_ = false;             //false: ion center and e- center overlap, true: there's a shift between the beam
    double bunch_separate_;
    int n_long_sample_ = 50;
public:
    IonSample ion_sample(){return ion_sample_;}
    unsigned int n_sample(){return n_sample_;}
    unsigned int n_tr(){return n_tr_;}
    unsigned int n_long(){return n_long_;}
    double bunch_separate(){return bunch_separate_;}
    bool shift(){return shift_;}
    int n_long_sample(){return n_long_sample_;}
    int set_shift(bool b){shift_ = b; return 0;}
    int set_n_tr(unsigned int n_tr){n_tr_ = n_tr; n_sample_=n_tr_*n_tr_*n_long_; return 0;}
    int set_n_long(unsigned int n){n_long_=n; n_sample_=n_tr_*n_tr_*n_long_; return 0;}
    int set_bunch_separate(double x){bunch_separate_ = x; return 0;}
    int set_n_long_sample(int x){ n_long_sample_ = x; return 0;}

    EcoolRateParas(){};
    EcoolRateParas(unsigned int n_sample):ion_sample_(IonSample::MONTE_CARLO),n_sample_(n_sample){};
    EcoolRateParas(unsigned int n_tr, unsigned int n_long):n_tr_(n_tr),n_long_(n_long){
        if (n_long_<2) n_long_ = 2;
        n_sample_ = n_tr_*n_tr_*n_long_;
        }
};
int ecooling_rate(EcoolRateParas &ecool_paras, ForceParas &force_paras, Beam &ion, Cooler &cooler, EBeam &ebeam,
                  Ring &ring, double &rate_x, double &rate_y, double &rate_s);
int end_ecooling(EcoolRateParas &ecool_paras, Beam &ion);
#endif // ECOOLING_H
