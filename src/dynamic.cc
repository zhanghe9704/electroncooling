#include "dynamic.h"
#include <chrono>
#include "constants.h"
#include "functions.h"
#include "particle_model.h"
#include "turn_by_turn.h"

DynamicParas *dynamic_paras = nullptr;
IBSParas *ibs_paras = nullptr;
EcoolRateParas *ecool_paras = nullptr;
ForceParas *force_paras = nullptr;

//std::unique_ptr<double []> rdn; //random number for model beam simulations

extern bool dynamic_flag;
extern double t_cooler;

extern std::unique_ptr<double []> x_bet, xp_bet, y_bet, yp_bet, ds, dp_p, x, y, xp, yp;
extern std::unique_ptr<double []> force_x, force_y, force_z;

int sample_the_ions(Beam &ion, Ring &ring, Cooler &cooler){
    switch (dynamic_paras->model()) {
    case DynamicModel::RMS : {
        if(dynamic_paras->ecool()) {
            config_ecooling(*ecool_paras, ion);
            ion_sample(*ecool_paras, ion, ring, cooler);
        }
        break;
    }
    case DynamicModel::TURN_BY_TURN : { //Same with PARTICLE model, except that coasting beam is sampled along the ring.
        initialize_turn_by_turn_model(ion, ring);
        break;
    }
    case DynamicModel::PARTICLE : {
        initialize_particle_model(ion);
        break;
    }
    default: {
        std::cout<<"Wrong dynamic model!"<<std::endl;
        assert(false);
    }
    }
    return 0;
}


int update_beam(int i, Beam &ion, Ring &ring, Cooler &cooler, EBeam &ebeam, std::vector<double> &r_ibs,
                std::vector<double> &r_ecool) {
    double dt = dynamic_paras->dt();
    double emit_nx = ion.emit_nx();
    double emit_ny = ion.emit_ny();
    double dp = ion.dp_p();
    double sigma_s = 0;
    if(ion.bunched()) sigma_s = ion.sigma_s();

    switch (dynamic_paras->model()) {
    case DynamicModel::RMS : {
        double rx = r_ibs.at(0) + r_ecool.at(0);
        double ry = r_ibs.at(1) + r_ecool.at(1);
        double rs = r_ibs.at(2) + r_ecool.at(2);

        //Calculate  new emittances
        emit_nx *= exp(rx*dt);
        emit_ny *= exp(ry*dt);
        dp *= dp*exp(rs*dt);
        dp = sqrt(dp);
        if(ion.bunched()) sigma_s = ring.beta_s()*dp;
        //update beam parameters
        ion.set_emit_nx(emit_nx);
        ion.set_emit_ny(emit_ny);
        ion.set_dp_p(dp);
        if(ion.bunched()) {
            if(dynamic_paras->fixed_bunch_length()) {
                ring.rf.v = ring.calc_rf_voltage();
            }
            else {
                ion.set_sigma_s(sigma_s);
            }
        }
//        if(ion.bunched()&& !dynamic_paras->fixed_bunch_length()) ion.set_sigma_s(sigma_s);
        //resample the ions
        if(dynamic_paras->ecool()) ion_sample(*ecool_paras, ion, ring, cooler);
        break;
    }
    case DynamicModel::PARTICLE : {
        if(dynamic_paras->ecool()) {
             double freq = k_c*ion.beta()/ring.circ()*cooler.section_number();
             //restore the coordinates
             restore_cord(t_cooler);
             //Adjust the frequency for bunched electron to coasting ion
             if(ebeam.bunched()&&(!ion.bunched())) {
                adjust_freq(freq, ebeam);
             }
            //Apply cooling kicks
            if(dynamic_paras->ecool()) {
                apply_cooling_kick(t_cooler, freq, dt);
            }
         }

        //Apply ibs kicks
        if(dynamic_paras->ibs()) {
            apply_ibs_kick(dt, ion, r_ibs);
        }

        move_particles(ion, ring);
//        //New betatron oscillation coordinates
//        double dx = dynamic_paras->twiss_ref.disp_x;
//        double dpx = dynamic_paras->twiss_ref.disp_dx;
//        double dy = dynamic_paras->twiss_ref.disp_y;
//        double dpy = dynamic_paras->twiss_ref.disp_dy;
//
//        for(unsigned int i=0; i<n_sample; ++i){
//            x_bet[i] = x[i] - dx*dp_p[i];
//            xp_bet[i] = xp[i] - dpx*dp_p[i];
//            y_bet[i] = y[i] - dy*dp_p[i];
//            yp_bet[i] = yp[i] - dpy*dp_p[i];
//        }
//        //random phase advance
//        double alf_x = dynamic_paras->twiss_ref.alf_x;
//        double alf_y = dynamic_paras->twiss_ref.alf_y;
//        double beta_x = dynamic_paras->twiss_ref.bet_x;
//        double beta_y = dynamic_paras->twiss_ref.bet_y;
//
//        double gamma_x = (1+alf_x*alf_x)/beta_x;
//        double gamma_y = (1+alf_y*alf_y)/beta_y;
//
//        uniform_random(n_sample, rdn.get(), -1, 1);
//        for(unsigned int i=0; i<n_sample; ++i){
//            double I = beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
//            double phi = k_pi*rdn[i];
//            x_bet[i] = sqrt(I*beta_x)*sin(phi);
//            xp_bet[i] = sqrt(I/beta_x)*(cos(phi)-alf_x*sin(phi));
//        }
//        uniform_random(n_sample, rdn.get(), -1, 1);
//        for(unsigned int i=0; i<n_sample; ++i){
//            double I = beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
//            double phi = k_pi*rdn[i];
//            y_bet[i] = sqrt(I*beta_y)*sin(phi);
//            yp_bet[i] = sqrt(I/beta_y)*(cos(phi)-alf_y*sin(phi));
//        }
//
//        if(ion.bunched()){
//            //temporarily use force_x to store the random numbers
//            uniform_random(n_sample, rdn.get(), -1, 1);
//            double beta_s = ring.beta_s();
//            double beta_s2_inv = 1/(beta_s*beta_s);
//            for(unsigned int i=0; i<n_sample; ++i){
//                double I = ds[i]*ds[i]*beta_s2_inv+dp_p[i]*dp_p[i];
//                I = sqrt(I);
//                double phi = k_pi*rdn[i];
//                dp_p[i] = I*sin(phi);
//                ds[i] = I*beta_s*cos(phi);
//            }
//        }
//
//        adjust_disp(dx, x_bet.get(), dp_p.get(), x.get(), n_sample);
//        adjust_disp(dy, y_bet.get(), dp_p.get(), y.get(), n_sample);
//        adjust_disp(dpx, xp_bet.get(), dp_p.get(), xp.get(), n_sample);
//        adjust_disp(dpy, yp_bet.get(), dp_p.get(), yp.get(), n_sample);

        //update beam parameters
        update_beam_parameters(ion);
        break;
    }
    case DynamicModel::TURN_BY_TURN : {
        //Apply ibs kicks
        if(dynamic_paras->ibs()) {
            apply_ibs_kick(dt, ion, r_ibs);
        }
        //Move particles
        turn_by_turn_move_particles(ion, ring, cooler);

        //update beam parameters
        update_beam_parameters(ion);
        break;
    }
    default: {
        assert(false&&"Wrong dynamic model!");
    }
    }
    return 0;
}

//void output(double t, std::vector<double> &emit, std::vector<double> &r, std::vector<double> &r_ibs,
//            std::vector<double> &r_ecool, bool bunched, std::ofstream &outfile) {
//    outfile<<t<<' '<<emit.at(0)<<' '<<emit.at(1)<<' '<<emit.at(2)<<' ';
//    if(bunched) outfile<<emit.at(3)<<' ';
//    else outfile<<0<<' ';
//    outfile<<r.at(0)<<' '<<r.at(1)<<' '<<r.at(2)<<' ';
//    outfile<<r_ibs.at(0)<<' '<<r_ibs.at(1)<<' '<<r_ibs.at(2)<<' ';
//    outfile<<r_ecool.at(0)<<' '<<r_ecool.at(1)<<' '<<r_ecool.at(2)<<' ';
//    outfile<<std::endl;
//}

void output(double t, std::vector<double> &emit, std::vector<double> &r, std::vector<double> &r_ibs,
            std::vector<double> &r_ecool, double v_rf, bool bunched, std::ofstream &outfile) {
    outfile<<t<<' '<<emit.at(0)<<' '<<emit.at(1)<<' '<<emit.at(2)<<' ';
    if(bunched) outfile<<emit.at(3)<<' ';
    else outfile<<0<<' ';
    outfile<<r.at(0)<<' '<<r.at(1)<<' '<<r.at(2)<<' ';
    outfile<<r_ibs.at(0)<<' '<<r_ibs.at(1)<<' '<<r_ibs.at(2)<<' ';
    outfile<<r_ecool.at(0)<<' '<<r_ecool.at(1)<<' '<<r_ecool.at(2)<<' ';
    outfile<<v_rf<<' ';
    outfile<<std::endl;
}

void save_ions_sdds(int n_sample, string filename) {
    using std::endl;
    std::ofstream output_particles;
    output_particles.open(filename);
    output_particles<<"SDDS1"<<endl;
    output_particles<<"! Define colums:"<<endl
        <<"&column name=x, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=xp, type=double, units=NULL, description=NULL, &end"<<endl
        <<"&column name=y, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=yp, type=double, units=NULL, description=NULL, &end"<<endl
        <<"&column name=ds, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=dp/p, type=double, units=NULL, description=NULL, &end"<<endl
        <<"!Declare ASCII data and end the header"<<endl
        <<"&data mode=ascii, &end"<<endl
        <<n_sample<<endl;
    output_particles.precision(10);
    output_particles<<std::showpos;
    output_particles<<std::scientific;
    for(int i=0; i<n_sample; ++i) {
        output_particles<<x[i]<<' '<<xp[i]<<' '<<y[i]<<' '<<yp[i]<<' '<<ds[i]<<' '<<dp_p[i]<<std::endl;
    }
    output_particles.close();
}


void save_ions(int n_sample, string filename) {
    std::ofstream output_particles;
    output_particles.open(filename);
    output_particles<<"x (m)              xp                 y (m)              yp                 ds (m)             dp_p"
        <<std::endl;
    output_particles.precision(10);
    output_particles<<std::showpos;
    output_particles<<std::scientific;
    for(int i=0; i<n_sample; ++i) {
        output_particles<<x[i]<<' '<<xp[i]<<' '<<y[i]<<' '<<yp[i]<<' '<<ds[i]<<' '<<dp_p[i]<<std::endl;
    }
    output_particles.close();
}


//void output_sddshead(int n, std::ofstream &outfile){
//    using std::endl;
//    outfile<<"SDDS1"<<endl;
//    outfile<<"! Define colums:"<<endl
//        <<"&column name=t, type=double, units=s, description=time, &end"<<endl
//        <<"&column name=emit_x, type=double, units=m*rad, description=\"normalized horizontal emittance\", &end"<<endl
//        <<"&column name=emit_y, type=double, units=m*rad, description=\"normalized vertical emittance\", &end"<<endl
//        <<"&column name=dp/p, type=double, units=NULL, description=\"momentum spread\", &end"<<endl
//        <<"&column name=sigma_s, type=double, units=m, description=\"RMS bunch length\", &end"<<endl
//        <<"&column name=rx, type=double, units=1/s, description=\"horizontal expansion rate\", &end"<<endl
//        <<"&column name=ry, type=double, units=1/s, description=\"vertical expansion rate\", &end"<<endl
//        <<"&column name=rs, type=double, units=1/s, description=\"longitudinal expansion rate\", &end"<<endl
//        <<"&column name=rx_ibs, type=double, units=1/s, description=\"horizontal IBS expansion rate\", &end"<<endl
//        <<"&column name=ry_ibs, type=double, units=1/s, description=\"vertical IBS expansion rate\", &end"<<endl
//        <<"&column name=rs_ibs, type=double, units=1/s, description=\"longitudinal IBS expansion rate\", &end"<<endl
//        <<"&column name=rx_ecool, type=double, units=1/s, description=\"horizontal electron cooling rate\", &end"<<endl
//        <<"&column name=ry_ecool, type=double, units=1/s, description=\"vertical electron cooling rate\", &end"<<endl
//        <<"&column name=rs_ecool, type=double, units=1/s, description=\"longitudinal electron cooling rate\", &end"<<endl
//        <<"!Declare ASCII data and end the header"<<endl
//        <<"&data mode=ascii, &end"<<endl
//        <<n<<endl;
//}

void output_sddshead(int n, std::ofstream &outfile){
    using std::endl;
    outfile<<"SDDS1"<<endl;
    outfile<<"! Define colums:"<<endl
        <<"&column name=t, type=double, units=s, description=time, &end"<<endl
        <<"&column name=emit_x, type=double, units=m*rad, description=\"normalized horizontal emittance\", &end"<<endl
        <<"&column name=emit_y, type=double, units=m*rad, description=\"normalized vertical emittance\", &end"<<endl
        <<"&column name=dp/p, type=double, units=NULL, description=\"momentum spread\", &end"<<endl
        <<"&column name=sigma_s, type=double, units=m, description=\"RMS bunch length\", &end"<<endl
        <<"&column name=rx, type=double, units=1/s, description=\"horizontal expansion rate\", &end"<<endl
        <<"&column name=ry, type=double, units=1/s, description=\"vertical expansion rate\", &end"<<endl
        <<"&column name=rs, type=double, units=1/s, description=\"longitudinal expansion rate\", &end"<<endl
        <<"&column name=rx_ibs, type=double, units=1/s, description=\"horizontal IBS expansion rate\", &end"<<endl
        <<"&column name=ry_ibs, type=double, units=1/s, description=\"vertical IBS expansion rate\", &end"<<endl
        <<"&column name=rs_ibs, type=double, units=1/s, description=\"longitudinal IBS expansion rate\", &end"<<endl
        <<"&column name=rx_ecool, type=double, units=1/s, description=\"horizontal electron cooling rate\", &end"<<endl
        <<"&column name=ry_ecool, type=double, units=1/s, description=\"vertical electron cooling rate\", &end"<<endl
        <<"&column name=rs_ecool, type=double, units=1/s, description=\"longitudinal electron cooling rate\", &end"<<endl
        <<"&column name=rf_voltage, type=double, units=V, description=\"Voltage of the RF cavity\", &end"<<endl
        <<"!Declare ASCII data and end the header"<<endl
        <<"&data mode=ascii, &end"<<endl
        <<n<<endl;
}

int dynamic(Beam &ion, Cooler &cooler, EBeam &ebeam, Ring &ring) {
    //Create temporary variables for expansion rate and beam macro-parameters
    std::vector<double> r_ibs = {0,0,0};
    std::vector<double> r_ecool = {0,0,0};
    std::vector<double> r = {0,0,0};
    std::vector<double> emit = {0,0,0,0} ; //emit_x, emit_y, dp_p, sigma_s

    bool ibs = dynamic_paras->ibs();
    bool ecool = dynamic_paras->ecool();
    int n_step = dynamic_paras->n_step();
    double dt = dynamic_paras->dt();
    double t = 0;
    int ion_save_itvl = dynamic_paras->ion_save_intvl();
    int output_itvl = dynamic_paras->output_intval();
    if (dynamic_paras->model()==DynamicModel::TURN_BY_TURN)
        dt = ring.circ()/(ion.beta()*k_c);

    //Prepare outputting
    std::ofstream outfile;
    outfile.open(dynamic_paras->output_file());
    output_sddshead(n_step+1, outfile);
    outfile.precision(10);
    outfile<<std::showpos;
    outfile<<std::scientific;

    //Set the twiss parameters for the reference point in particle model simulation
    referece_twiss(cooler);

    if(!ecool && dynamic_paras->model()==DynamicModel::PARTICLE)
        assert(dynamic_paras->n_sample()>0 && "Sample particle number must be greater than zero!");

    //sample the ion
    sample_the_ions(ion, ring, cooler);

//    if(dynamic_paras->fixed_bunch_length()) {
//        ring.rf.v = ring.calc_rf_voltage();
//    }

    if(ring.rf.gamma_tr>0) ring.rf.v = ring.calc_rf_voltage();

    //Start tracking
    std::cout<<"Start dynamic simulation ... "<<std::endl;
    dynamic_flag = true;

    for(int i=0; i<n_step+1; ++i) {
        if (ion_save_itvl>0 && i%ion_save_itvl==0 && dynamic_paras->model()!=DynamicModel::RMS)
            save_ions_sdds(dynamic_paras->n_sample(), "ions"+std::to_string(i)+".txt");

        //record
        emit.at(0) = ion.emit_nx();
        emit.at(1) = ion.emit_ny();
        emit.at(2) = ion.dp_p();
        if (ion.bunched()) emit.at(3) = ion.sigma_s();

        //Rate calculation
        if(ibs) ibs_rate(*ring.lattice_, ion, *ibs_paras, r_ibs.at(0), r_ibs.at(1), r_ibs.at(2));
        if(ecool) {
            ecooling_rate(*ecool_paras, *force_paras, ion, cooler, ebeam, ring, r_ecool.at(0), r_ecool.at(1), r_ecool.at(2));
        }
        for(int i=0; i<3; ++i) r.at(i) = r_ibs.at(i) + r_ecool.at(i);

        //Output
        if (output_itvl==1) output(t, emit, r, r_ibs, r_ecool, ring.rf.v, ion.bunched(), outfile);
        else if(i%output_itvl==0) output(t, emit, r, r_ibs, r_ecool, ring.rf.v, ion.bunched(), outfile);

        //Update beam parameters and particles
        update_beam(i, ion, ring, cooler, ebeam, r_ibs, r_ecool);
        t += dt;
        std::cout<<i<<std::endl;
    }

    if (ion_save_itvl>0 && n_step%ion_save_itvl!=0 && dynamic_paras->model()!=DynamicModel::RMS)
        save_ions_sdds(dynamic_paras->n_sample(), "ions"+std::to_string(n_step)+".txt");
    dynamic_flag = false;
    outfile.close();
    std::cout<<"Finished dynamic simulation."<<std::endl;

    return 0;
}


