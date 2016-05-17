#ifndef IBS_HPP
#define IBS_HPP

#include "ring.h"
#include <cmath>
#include <cstring>
#include <iostream>

class IBSParas{
    int nu_ = 0;                //Grid number in u direction.
    int nv_ = 0;                //Grid number in v direction.
    int nz_ = 0;                //Grid number in z direction.
    double log_c_ = 20;     //Coulomb logarithm.
    double k_ = 0;          //Coupling rate in transverse directions.
    bool use_log_c_ = true;
    bool reset_ = true;
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
    IBSParas(int nu, int nv):nu_(nu),nv_(nv){};
    IBSParas(int nu, int nv, double log_c):nu_(nu),nv_(nv),log_c_(log_c){};
    IBSParas(int nu, int nv, int nz):nu_(nu),nv_(nv),nz_(nz){use_log_c_ = false;}

};

int config_ibs(Lattice &lattice);
int ibs_rate(Lattice &lattice, Beam &beam, IBSParas &ibs_paras,double &rx, double &ry, double &rs);
int end_ibs();

#endif
