#include "turn_by_turn.h"
#include <chrono>
#include <cmath>
#include "beam.h"
#include "constants.h"
#include "dynamic.h"
#include "ecooling.h"
#include "functions.h"
#include "particle_model.h"
#include "ring.h"

extern DynamicParas *dynamic_paras;
extern EcoolRateParas *ecool_paras;
extern std::unique_ptr<double []> x_bet, xp_bet, y_bet, yp_bet, ds, dp_p, x, y, xp, yp;
extern std::unique_ptr<double []> force_x, force_y, force_z;

extern std::unique_ptr<double []> rdn;
extern unsigned int n_sample;
extern double p0;

void initialize_turn_by_turn_model(Beam &ion, Ring &ring) {
    n_sample = dynamic_paras->n_sample();
    p0 = ion.gamma()*ion.mass()*1e6*k_e*ion.beta()/k_c;
    if(dynamic_paras->ecool()) {
        assert(ecool_paras->ion_sample()==IonSample::MONTE_CARLO);
        config_ecooling(*ecool_paras, ion);
    }
    else {
        assert(n_sample>0&&"The number of the ions should be greater than zero! Define n_sample in Dynamic_paras.");
        particle_model_ibs_scratches();
    }
    assert(dynamic_paras->twiss_ref.bet_x>0&& dynamic_paras->twiss_ref.bet_y>0
           &&"Need to define the TWISS parameters for the reference point for simulations using the turn-by-turn model!");
    ion_beam_model_MonteCarlo_Gaussian(n_sample, ion, dynamic_paras->twiss_ref);

    if(!ion.bunched()) {      //Sample the coating ion beam around the ring
        uniform_random(n_sample, ds.get(), -1*ring.circ()/2, ring.circ()/2);
        uniform_random_adjust(n_sample, ds.get());
        gaussian_random(n_sample, dp_p.get(), ion.dp_p());
        gaussian_random_adjust(n_sample, dp_p.get(),ion.dp_p());
    }
    if (rdn.get() == nullptr) rdn.reset(new double[n_sample]);
}
////
void turn_by_turn_move_particles(Beam &ion, Ring &ring, Cooler &cooler) {
    //Transverse
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
    //one turn phase advance
//    double alf_x = dynamic_paras->twiss_ref.alf_x;
//    double alf_y = dynamic_paras->twiss_ref.alf_y;
//    double beta_x = dynamic_paras->twiss_ref.bet_x;
//    double beta_y = dynamic_paras->twiss_ref.bet_y;

//    double gamma_x = (1+alf_x*alf_x)/beta_x;
//    double gamma_y = (1+alf_y*alf_y)/beta_y;

    //Transverse motion by tunes
    assert(ring.tunes.qx>0&&ring.tunes.qy>0&&"Transverse tunes are needed for Turn_by_turn model");
    double Qx = ring.tunes.qx;
    double Qy = ring.tunes.qy;
    for (unsigned int i=0; i<n_sample; ++i) {
        double phi = 2*k_pi*Qx;
        double x_bet_0 = x_bet[i];
        double xp_bet_0 = xp_bet[i];
        x_bet[i] = cos(phi)*x_bet_0 + cooler.beta_h()*sin(phi)*xp_bet_0;
        xp_bet[i] = -sin(phi)*x_bet_0/cooler.beta_h() + cos(phi)*xp_bet_0;
        phi = 2*k_pi*Qy;
        double y_bet_0 = y_bet[i];
        double yp_bet_0 = yp_bet[i];
        y_bet[i] = cos(phi)*y_bet_0 + sin(phi)*yp_bet_0*cooler.beta_v();
        yp_bet[i] = -sin(phi)*y_bet_0/cooler.beta_v() + cos(phi)*yp_bet_0;
    }

    //Longitudinal motion.
    if (ring.tunes.qs>0||ring.rf.v>0) {    //RF, synchrotron oscillation.
//        assert(ring.tunes->qs>0||ring.rf->v>0&&"Longitudinal tune or RF cavity needed for Turn_by_turn model");

        if(ring.rf.v>0) { //Longitudinal motion by RF.
            double circ = ring.circ();
            double beta2 = ion.beta()*ion.beta();
            double beta2_inv = 1/beta2;
            double total_energy = ion.gamma()*ion.mass(); //ion total energy [MeV/c^2]
            double total_energy_inv = 1/total_energy;
            double adj_dp2dE = beta2*total_energy;

            double volt = ring.rf.v;
            double phi_s = ring.rf.phi;
//            double phi_0 = ring.rf->phi_0();
            double h = ring.rf.h;
//            double s_s = phi_s*circ/(h*2*k_pi);
            double half_phase = h*k_pi;
            double total_phase = h*2*k_pi;
            double adj_s2phi = total_phase/circ;
            double adj_phi2s = 1/adj_s2phi;
            double adj_dE = ion.charge_number()*volt*1e-6; // [MeV/c^2]
    //                double adj_dE = ion.charge_number()*ring.rf_->volt()*1e-6; // [MeV/c^2]
            double sin_phi_s = sin(phi_s);
//            double sin_phi_s = sin(phi_s+phi_0);
            double eta = 1/(ring.rf.gamma_tr*ring.rf.gamma_tr) - 1/(ion.gamma()*ion.gamma()); //phase slip factor
            double adj_dE2dphi = total_phase*eta*beta2_inv*total_energy_inv;
            for(unsigned int i = 0; i < n_sample; ++i) {
                dp_p[i] *= adj_dp2dE; //dp/p -> dE/E -> dE in [MeV/c^2]
//                ds[i] += s_s;  //s = ds + s_s: adjust ds to be measured from the start of the ring
                ds[i] *= adj_s2phi;  //phi = s*h*2*pi/circ: s -> phi
//                dp_p[i] += adj_dE*(sin(ds[i]+phi_0)-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                dp_p[i] += adj_dE*(sin(ds[i])-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                ds[i] += adj_dE2dphi*dp_p[i]; //phi_n+1 = phi_n + 2*pi*h*eta*dE_n+1 / (beta^2 * E)
                ds[i] += half_phase;    //phi in [0, total phase]
                ds[i] = fmod(ds[i], total_phase);  //phi_n+1 module h*2*pi
                if (ds[i]<0) ds[i] += total_phase; //adjust phi is phi is less than zero
                ds[i] -= half_phase;    //phi in [-half_phase, half_phase]
                ds[i] *= adj_phi2s;  //phi -> s
//                ds[i] -= s_s;        // ds = s - s_s: adjust ds back to be centered about s_s
                dp_p[i]  = dp_p[i]*total_energy_inv*beta2_inv; //dE -> dE/E -> dp/p = beta*beta*dE/E;
            }
        }
        else if(ring.tunes.qs>0) {//Longitudinal motion by tune
            double phi = 2*k_pi*ring.tunes.qs;
            double inv_beta_s = 1/ring.beta_s();
            double beta_s = ring.beta_s();
            for (unsigned int i=0; i<n_sample; ++i) {
                double dp_p_0 = dp_p[i];
                double ds_0 = ds[i];
                dp_p[i] = cos(phi)*dp_p_0 - sin(phi)*ds_0*inv_beta_s;
                ds[i] = sin(phi)*dp_p_0*beta_s + cos(phi)*ds_0;
            }
        }
    }
    else {  //No RF.
        double gamma_0 = ion.gamma();
        double beta_0 = ion.beta();
        double half_length = 0.5*ring.circ();
        for(unsigned int i=0; i<n_sample; ++i) {
            double gamma2 = 1+(1+dp_p[i])*(1+dp_p[i])*(gamma_0*gamma_0-1);
            double beta = sqrt(1-1/gamma2);
            double s = (beta/beta_0-1)*2*half_length;
            ds[i] += half_length;    //s in [0, 2*half_length]
            ds[i] += s;
            ds[i] = fmod(ds[i], 2*half_length);
            if(ds[i]<0) ds[i] += 2*half_length;
            ds[i] -= half_length;   //s in [-half_length, half_length]
        }
    }

    adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
    adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
    adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
    adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);

}
