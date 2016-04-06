#include "dynamic.h"

DynamicParas *dynamic_paras = nullptr;
IBSParas *ibs_paras = nullptr;
EcoolRateParas *ecool_paras = nullptr;
ForceParas *force_paras = nullptr;

double *rcd_emit_x;
double *rcd_emit_y;
double *rcd_sigma_s;
double *rcd_dp_p;
double *rcd_ibs_rate_x, *rcd_ibs_rate_y, *rcd_ibs_rate_s, *rcd_ecool_rate_x, *rcd_ecool_rate_y, *rcd_ecool_rate_s;

extern int model_beam_count;
extern int rms_dynamic_count;
extern double *x_bet, *xp_bet, *y_bet, *yp_bet, *ds, *dp_p, *x, *y, *xp, *yp;

int record_config(int n, bool bunched) {
    int n_step = n+1;
    rcd_emit_x = new double[n_step];
    rcd_emit_y = new double[n_step];
    rcd_dp_p = new double[n_step];
    rcd_ibs_rate_x = new double[n_step];
    rcd_ibs_rate_y = new double[n_step];
    rcd_ibs_rate_s = new double[n_step];
    rcd_ecool_rate_x = new double[n_step];
    rcd_ecool_rate_y = new double[n_step];
    rcd_ecool_rate_s = new double[n_step];
    if (bunched) rcd_sigma_s = new double[n_step];
    return 0;
}

int clear_record(bool bunched) {
    delete[] rcd_emit_x;
    delete[] rcd_emit_y;
    delete[] rcd_dp_p;
    delete[] rcd_ibs_rate_x;
    delete[] rcd_ibs_rate_y;
    delete[] rcd_ibs_rate_s;
    delete[] rcd_ecool_rate_x;
    delete[] rcd_ecool_rate_y;
    delete[] rcd_ecool_rate_s;
    if (bunched) delete[] rcd_sigma_s;
    return 0;
}

int dynamic_backup(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring,std::ofstream &outfile) {
    double rx_ibs, ry_ibs, rs_ibs, rx_ecool, ry_ecool, rs_ecool, rx, ry, rs;
    bool ibs = dynamic_paras->ibs();
    bool ecool = dynamic_paras->ecool();
    double bs = 0;
    if (ion.bunched()) bs = ring.beta_s();

    rx_ibs = 0;
    ry_ibs = 0;
    rs_ibs = 0;
    rx_ecool = 0;
    ry_ecool = 0;
    rs_ecool = 0;
    double dt = dynamic_paras->dt();
    double emit_nx = ion.emit_nx();
    double emit_ny = ion.emit_ny();
    double dp = ion.dp_p();
    double sigma_s = 0;
    if(ion.bunched()) sigma_s = ion.sigma_s();
    double t = 0;
    if(ibs) config_ibs(*ring.lattice_);
    for(int i=0; i<dynamic_paras->n_step()+1; ++i) {
        if(ibs) ibs_rate(*ring.lattice_, ion, *ibs_paras, rx_ibs, ry_ibs, rs_ibs);
        if(ecool) {
            ++rms_dynamic_count;
            ecooling_rate(*ecool_paras, *force_paras, ion, cooler, ebeam, ring, rx_ecool, ry_ecool, rs_ecool);
        }
        rx = rx_ibs + rx_ecool;
        ry = ry_ibs + ry_ecool;
        rs = rs_ibs + rs_ecool;

        emit_nx *= exp(rx*dt);
        emit_ny *= exp(ry*dt);
        dp *= dp*exp(rs*dt);
        dp = sqrt(dp);
        if(ion.bunched()) sigma_s = bs*dp;
        t += dt;

        ion.set_emit_nx(emit_nx);
        ion.set_emit_ny(emit_ny);
        ion.set_dp_p(dp);
        if(ion.bunched()) ion.set_sigma_s(sigma_s);

        //output to file
    }

    if(ibs) end_ibs();
    if(ecool) {
        rms_dynamic_count = -1;
        end_ecooling(*ecool_paras, ion);
    }
    return 0;

}

int dynamic(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring, std::ofstream &outfile) {
    //Initialization

    double rx_ibs, ry_ibs, rs_ibs, rx_ecool, ry_ecool, rs_ecool, rx, ry, rs;
    bool ibs = dynamic_paras->ibs();
    bool ecool = dynamic_paras->ecool();
    double bs = 0;
    if (ion.bunched()) bs = ring.beta_s();

    rx_ibs = 0;
    ry_ibs = 0;
    rs_ibs = 0;
    rx_ecool = 0;
    ry_ecool = 0;
    rs_ecool = 0;
    int n_step = dynamic_paras->n_step();
    double dt = dynamic_paras->dt();
    double emit_nx = ion.emit_nx();
    double emit_ny = ion.emit_ny();
    double dp = ion.dp_p();
    double sigma_s = 0;
    if(ion.bunched()) sigma_s = ion.sigma_s();
    double t = 0;
    if(ibs) config_ibs(*ring.lattice_);

    record_config(n_step, ion.bunched());

    //sample the ion
    //For now the sample is done in the first run of the cooling rate calculation.
    //Probably it should be move out of the loop to here.

    //Start tracking
    std::cout<<"Start dynamic simulation ... "<<std::endl;
    for(int i=0; i<n_step+1; ++i) {
        //record
        rcd_emit_x[i] = emit_nx;
        rcd_emit_y[i] = emit_ny;
        rcd_dp_p[i] = dp;
        if (ion.bunched()) rcd_sigma_s[i] = sigma_s;


        //IBS rate
        if(ibs) ibs_rate(*ring.lattice_, ion, *ibs_paras, rx_ibs, ry_ibs, rs_ibs);
        rcd_ibs_rate_x[i] = rx_ibs;
        rcd_ibs_rate_y[i] = ry_ibs;
        rcd_ibs_rate_s[i] = rs_ibs;
        //Cooling rate
        if(ecool) {
            ++rms_dynamic_count;
            if(dynamic_paras->model()==DynamicModel::MODEL_BEAM) ++model_beam_count;
            ecooling_rate(*ecool_paras, *force_paras, ion, cooler, ebeam, ring, rx_ecool, ry_ecool, rs_ecool);

        }
        rcd_ecool_rate_x[i] = rx_ecool;
        rcd_ecool_rate_y[i] = ry_ecool;
        rcd_ecool_rate_s[i] = rs_ecool;

        rx = rx_ibs + rx_ecool;
        ry = ry_ibs + ry_ecool;
        rs = rs_ibs + rs_ecool;

        //Output
        outfile<<t<<' '<<rcd_emit_x[i]<<' '<<rcd_emit_y[i]<<' '<<rcd_dp_p[i]<<' ';
        if(ion.bunched()) outfile<<rcd_sigma_s[i]<<' ';
        else outfile<<0<<' ';
        outfile<<rx<<' '<<ry<<' '<<rs<<' ';
        outfile<<std::endl;

        //Update particles
        switch (dynamic_paras->model()) {
        case DynamicModel::RMS:
            emit_nx *= exp(rx*dt);
            emit_ny *= exp(ry*dt);
            dp *= dp*exp(rs*dt);
            dp = sqrt(dp);
            if(ion.bunched()) sigma_s = bs*dp;
            break;
        case DynamicModel::MODEL_BEAM:
            break;
        }
        t += dt;

        ion.set_emit_nx(emit_nx);
        ion.set_emit_ny(emit_ny);
        ion.set_dp_p(dp);
        if(ion.bunched()) ion.set_sigma_s(sigma_s);

        std::cout<<i<<std::endl;
    }
    std::cout<<"Finished dynamic simulation."<<std::endl;

    //Clean;
    clear_record(ion.bunched());
    if(ibs) end_ibs();
    if(ecool) {
        rms_dynamic_count = -1;
        end_ecooling(*ecool_paras, ion);
    }
    return 0;

}


