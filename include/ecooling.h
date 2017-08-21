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
    int n_sample_ = 125000;
    int n_tr_ = 50;
    int n_long_ = 50;
    bool shift_ = false;             //false: ion center and e- center overlap, true: there's a shift between the beam
    double bunch_separate_;
    int n_long_sample_ = 50;
public:
    IonSample ion_sample(){return ion_sample_;}
    int n_sample(){return n_sample_;}
    int n_tr(){return n_tr_;}
    int n_long(){return n_long_;}
    double bunch_separate(){return bunch_separate_;}
    bool shift(){return shift_;}
    int n_long_sample(){return n_long_sample_;}
    int set_n_sample(int n_sample){n_sample_ = n_sample; return 0;}
    int set_shift(bool b){shift_ = b; return 0;}
    int set_n_tr(unsigned int n_tr){n_tr_ = n_tr; n_sample_=n_tr_*n_tr_*n_long_; return 0;}
    int set_n_long(unsigned int n){n_long_=n; n_sample_=n_tr_*n_tr_*n_long_; return 0;}
    int set_bunch_separate(double x){bunch_separate_ = x; return 0;}
    int set_n_long_sample(int x){ n_long_sample_ = x; return 0;}

    EcoolRateParas(){};
    EcoolRateParas(int n_sample):ion_sample_(IonSample::MONTE_CARLO),n_sample_(n_sample){};
    EcoolRateParas(int n_tr, int n_long):n_tr_(n_tr),n_long_(n_long){
        if (n_long_<2) n_long_ = 2;
        n_sample_ = n_tr_*n_tr_*n_long_;
        }
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
//    Twiss():bet_x_(0),bet_y_(0),alf_x_(0),alf_y_(0),disp_x_(0),disp_y_(0),disp_dx_(0),disp_dy_(0){};
};

int ecooling_rate(EcoolRateParas &ecool_paras, ForceParas &force_paras, Beam &ion, Cooler &cooler, EBeam &ebeam,
                  Ring &ring, double &rate_x, double &rate_y, double &rate_s);
//int end_ecooling(EcoolRateParas &ecool_paras, Beam &ion);

int config_ecooling(EcoolRateParas &ecool_paras, Beam &ion);
int ion_sample(EcoolRateParas &ecool_paras, Beam &ion, Ring &ring, Cooler &cooler);
int ion_beam_model_MonteCarlo_Gaussian(unsigned int n_sample, Beam &ion, Twiss &twiss);
double emit(double * x, double * xp, unsigned int n);
double emit_p(double * dp_p, unsigned int n);
double emit_p(double * dp_p, double * ds, Ring &ring, unsigned int n);
int adjust_disp(double dx, double *x_bet, double *dp_p, double *x, unsigned int n);
#endif // ECOOLING_H
