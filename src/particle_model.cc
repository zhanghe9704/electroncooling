#include "particle_model.h"
#include <chrono>
#include <cmath>
#include "constants.h"
#include "dynamic.h"
#include "ecooling.h"
#include "functions.h"

extern DynamicParas *dynamic_paras;
extern EcoolRateParas *ecool_paras;
extern std::unique_ptr<double []> x_bet, xp_bet, y_bet, yp_bet, ds, dp_p, x, y, xp, yp;
extern std::unique_ptr<double []> force_x, force_y, force_z;

std::unique_ptr<double []> rdn;
unsigned int n_sample = 0;
double p0 = 0;

void initialize_particle_model(Beam &ion) {
    n_sample = dynamic_paras->n_sample();
    p0 = ion.gamma()*ion.mass()*1e6*k_e*ion.beta()/k_c;
    if(dynamic_paras->ecool()) {
        assert(ecool_paras->ion_sample()==IonSample::MONTE_CARLO);
        config_ecooling(*ecool_paras, ion);
//        n_sample = ecool_paras->n_sample();
    }
    else {
        assert(n_sample>0&&"The number of the ions should be greater than zero! Define n_sample in Dynamic_paras.");
        particle_model_ibs_scratches();
    }
    assert(dynamic_paras->twiss_ref.bet_x>0&& dynamic_paras->twiss_ref.bet_y>0
           &&"Need to define the TWISS parameters for the reference point for simulations using the particle model!");
    ion_beam_model_MonteCarlo_Gaussian(n_sample, ion, dynamic_paras->twiss_ref);
    if (rdn.get() == nullptr) rdn.reset(new double[n_sample]);
}

void referece_twiss(Cooler &cooler) {
    if(dynamic_paras->model()==DynamicModel::PARTICLE||dynamic_paras->model()==DynamicModel::TURN_BY_TURN) {
        if(dynamic_paras->ecool()) {
            //set the reference point at the cooler
            dynamic_paras->twiss_ref.bet_x = cooler.beta_h();
            dynamic_paras->twiss_ref.bet_y = cooler.beta_v();
            dynamic_paras->twiss_ref.alf_x = cooler.alpha_h();
            dynamic_paras->twiss_ref.alf_y = cooler.alpha_v();
            dynamic_paras->twiss_ref.disp_x = cooler.disp_h();
            dynamic_paras->twiss_ref.disp_y = cooler.disp_v();
            dynamic_paras->twiss_ref.disp_dx = cooler.der_disp_h();
            dynamic_paras->twiss_ref.disp_dy = cooler.der_disp_v();
        }
    }
}

//int model_beam_ibs_scratches(int n_sample) {
void particle_model_ibs_scratches() {
    x_bet.reset(new double[n_sample]);
    xp_bet.reset(new double[n_sample]);
    y_bet.reset(new double[n_sample]);
    yp_bet.reset(new double[n_sample]);
    ds.reset(new double[n_sample]);
    dp_p.reset(new double[n_sample]);
    x.reset(new double[n_sample]);
    y.reset(new double[n_sample]);
    xp.reset(new double[n_sample]);
    yp.reset(new double[n_sample]);
    srand(time(NULL));
}

void ibs_kick(unsigned int n_sample, double rate, double twiss, double dt, double emit, double* p) {

    if (rate>0) {
        double theta = sqrt(2*rate*dt*emit/twiss);
        gaussian_random(n_sample, rdn.get(), 1, 0);
        for(unsigned int i=0; i<n_sample; ++i) p[i] += theta*rdn[i];
    }
    else {
        double k = exp(rate*dt);
        for(unsigned int i=0; i<n_sample; ++i) p[i] *= k;
    }
}

void restore_cord(double t_cooler) {
    double inv_p0 = 1.0/p0;
    for(unsigned int i=0; i<n_sample; ++i) {
         xp[i] -= force_x[i]*t_cooler*inv_p0;
         yp[i] -= force_y[i]*t_cooler*inv_p0;
         dp_p[i] -= force_z[i]*t_cooler*inv_p0;
     }
}

void adjust_freq(double &freq, EBeam ebeam) {
    double sample_length = ebeam.length();
    double bunch_separate = ecool_paras->bunch_separate();
    if(sample_length<0) perror("electron bunch length must be positive!");
    if(bunch_separate>sample_length) {
        double duty_factor = sample_length/bunch_separate;
        freq *= duty_factor;
    }
    else {
        perror("Electron bunch length is larger than the distance between electron bunches");
    }
}

void apply_cooling_kick(double t_cooler, double freq, double dt) {
    for(unsigned int i=0; i<n_sample; ++i) {
        xp[i] = !iszero(xp[i])?xp[i]*exp(force_x[i]*t_cooler*dt*freq/(xp[i]*p0)):xp[i];
        yp[i] = !iszero(yp[i])?yp[i]*exp(force_y[i]*t_cooler*dt*freq/(yp[i]*p0)):yp[i];
        dp_p[i] = !iszero(dp_p[i])?dp_p[i]*exp(force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0)):dp_p[i];
//        xp[i] = xp[i]!=0?xp[i]*exp(force_x[i]*t_cooler*dt*freq/(xp[i]*p0)):xp[i];
//        yp[i] = yp[i]!=0?yp[i]*exp(force_y[i]*t_cooler*dt*freq/(yp[i]*p0)):yp[i];
//        dp_p[i] = dp_p[i]!=0?dp_p[i]*exp(force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0)):dp_p[i];
    }
}

void apply_ibs_kick(double dt, Beam &ion, std::vector<double> &r_ibs) {
    assert(dynamic_paras->twiss_ref.bet_x>0&& dynamic_paras->twiss_ref.bet_y>0
                   &&"TWISS parameters for the reference point not defined! Define twiss_ref.");
    double rx_ibs = r_ibs.at(0);
    double ry_ibs = r_ibs.at(1);
    double rs_ibs = r_ibs.at(2);

    ibs_kick(n_sample, rx_ibs, dynamic_paras->twiss_ref.bet_x, dt, ion.emit_x(), xp.get());
    ibs_kick(n_sample, ry_ibs, dynamic_paras->twiss_ref.bet_y, dt, ion.emit_y(), yp.get());
    if (ion.bunched())
        ibs_kick(n_sample, rs_ibs, 1, dt, ion.dp_p()*ion.dp_p(), dp_p.get());
    else
        ibs_kick(n_sample, rs_ibs, 2, dt, ion.dp_p()*ion.dp_p(), dp_p.get());
}

void move_particles(Beam &ion, Ring &ring) {
    //New betatron oscillation coordinates
    double dx = dynamic_paras->twiss_ref.disp_x;
    double dpx = dynamic_paras->twiss_ref.disp_dx;
    double dy = dynamic_paras->twiss_ref.disp_y;
    double dpy = dynamic_paras->twiss_ref.disp_dy;

    for(unsigned int i=0; i<n_sample; ++i){
        x_bet[i] = x[i] - dx*dp_p[i];
        xp_bet[i] = xp[i] - dpx*dp_p[i];
        y_bet[i] = y[i] - dy*dp_p[i];
        yp_bet[i] = yp[i] - dpy*dp_p[i];
    }
    //random phase advance
    double alf_x = dynamic_paras->twiss_ref.alf_x;
    double alf_y = dynamic_paras->twiss_ref.alf_y;
    double beta_x = dynamic_paras->twiss_ref.bet_x;
    double beta_y = dynamic_paras->twiss_ref.bet_y;

    double gamma_x = (1+alf_x*alf_x)/beta_x;
    double gamma_y = (1+alf_y*alf_y)/beta_y;

    uniform_random(n_sample, rdn.get(), -1, 1);
    for(unsigned int i=0; i<n_sample; ++i){
        double I = beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
        double phi = k_pi*rdn[i];
        x_bet[i] = sqrt(I*beta_x)*sin(phi);
        xp_bet[i] = sqrt(I/beta_x)*(cos(phi)-alf_x*sin(phi));
    }
    uniform_random(n_sample, rdn.get(), -1, 1);
    for(unsigned int i=0; i<n_sample; ++i){
        double I = beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
        double phi = k_pi*rdn[i];
        y_bet[i] = sqrt(I*beta_y)*sin(phi);
        yp_bet[i] = sqrt(I/beta_y)*(cos(phi)-alf_y*sin(phi));
    }

    if(ion.bunched()){
        //temporarily use force_x to store the random numbers
        uniform_random(n_sample, rdn.get(), -1, 1);
        double beta_s = ring.beta_s();
        double beta_s2_inv = 1/(beta_s*beta_s);
        for(unsigned int i=0; i<n_sample; ++i){
            double I = ds[i]*ds[i]*beta_s2_inv+dp_p[i]*dp_p[i];
            I = sqrt(I);
            double phi = k_pi*rdn[i];
            dp_p[i] = I*sin(phi);
            ds[i] = I*beta_s*cos(phi);
        }
    }

    adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
    adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
    adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
    adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);
}

void update_beam_parameters(Beam &ion) {
    double emit_x = emit(x_bet.get(), xp_bet.get(), n_sample);
    double emit_y = emit(y_bet.get(), yp_bet.get(), n_sample);
    double dp = sqrt(emit_p(dp_p.get(), n_sample));

    ion.set_emit_x(emit_x);
    ion.set_emit_y(emit_y);
    ion.set_dp_p(dp);

    if(ion.bunched()) {
        double sigma_s = sqrt(emit_p(ds.get(), n_sample));
        ion.set_sigma_s(sigma_s);
    }
}
