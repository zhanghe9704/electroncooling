#ifndef IBS_HPP
#define IBS_HPP

#include "ring.h"
#include <assert.h>
#include <cmath>
#include <cstring>
#include <iostream>

enum class IBSModel {MARTINI, BM};

class IBSParas{
    int nu_ = 0;                //Grid number in u direction.
    int nv_ = 0;                //Grid number in v direction.
    int nz_ = 0;                //Grid number in z direction.
    double log_c_ = 20;     //Coulomb logarithm.
    double k_ = 0;          //Coupling rate in transverse directions.
    bool use_log_c_ = true;
    bool reset_ = true;
    IBSModel model_ = IBSModel::MARTINI;
public:
    int nu(){return nu_;}
    int nv(){return nv_;}
    int nz(){return nz_;}
    double log_c(){return log_c_;}
    double k(){return k_;}
    bool use_log_c(){return use_log_c_;}
    bool reset(){return reset_;}
    int set_k(double x){k_ = x; return 0;}
    int set_log_c(double x){ log_c_ = x; use_log_c_ = true; return 0;}
    int reset_off(){reset_ = false; return 0;}
    int reset_on(){reset_ = true; return 0;}
    void set_model(IBSModel model) {model_ = model;}
    void set_nu(int nu){assert(nu>0&&"Wrong value of nu in IBS parameters!"); nu_ = nu;}
    void set_nv(int nv){assert(nv>0&&"Wrong value of nv in IBS parameters!"); nv_ = nv;}
    void set_nz(int nz){assert(nz>0&&"Wrong value of nz in IBS parameters!"); nz_ = nz;}
    IBSModel model() {return model_;}

    IBSParas(int nu, int nv):nu_(nu),nv_(nv) { }
    IBSParas(int nu, int nv, double log_c):nu_(nu),nv_(nv),log_c_(log_c) { }
    IBSParas(int nu, int nv, int nz):nu_(nu),nv_(nv),nz_(nz) { use_log_c_ = false; }
    IBSParas(IBSModel model):model_(model) { }
};

void ibs_rate(Lattice &lattice, Beam &beam, IBSParas &ibs_paras,double &rx, double &ry, double &rs);


#endif
