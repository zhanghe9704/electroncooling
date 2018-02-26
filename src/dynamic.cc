#include "dynamic.h"
#include <chrono>
#include "constants.h"
#include "particle_model.h"
#include "functions.h"

DynamicParas *dynamic_paras = nullptr;
IBSParas *ibs_paras = nullptr;
EcoolRateParas *ecool_paras = nullptr;
ForceParas *force_paras = nullptr;

//std::unique_ptr<double []> rdn; //random number for model beam simulations

extern bool dynamic_flag;
extern double t_cooler;

extern std::unique_ptr<double []> x_bet, xp_bet, y_bet, yp_bet, ds, dp_p, x, y, xp, yp;
extern std::unique_ptr<double []> force_x, force_y, force_z;

////int model_beam_ibs_scratches(int n_sample) {
//int particle_model_ibs_scratches(int n_sample) {
//    x_bet.reset(new double[n_sample]);
//    xp_bet.reset(new double[n_sample]);
//    y_bet.reset(new double[n_sample]);
//    yp_bet.reset(new double[n_sample]);
//    ds.reset(new double[n_sample]);
//    dp_p.reset(new double[n_sample]);
//    x.reset(new double[n_sample]);
//    y.reset(new double[n_sample]);
//    xp.reset(new double[n_sample]);
//    yp.reset(new double[n_sample]);
//    srand(time(NULL));
//    return 0;
//}

int sample_the_ions(Beam &ion, Ring &ring, Cooler &cooler){
    switch (dynamic_paras->model()) {
    case DynamicModel::RMS : {
        if(dynamic_paras->ecool()) {
            config_ecooling(*ecool_paras, ion);
            ion_sample(*ecool_paras, ion, ring, cooler);
        }
        break;
    }
//    case DynamicModel::MODEL_BEAM : {}
    case DynamicModel::PARTICLE : {
        initialize_particle_model(ion);
//        int n_sample = dynamic_paras->n_sample();
//        if(dynamic_paras->ecool()) {
//            assert(ecool_paras->ion_sample()==IonSample::MONTE_CARLO);
//            config_ecooling(*ecool_paras, ion);
//            n_sample = ecool_paras->n_sample();
//        }
//        else {
//            assert(n_sample>0&&"The number of the ions should be greater than zero! Define n_sample in Dynamic_paras.");
//            particle_model_ibs_scratches(n_sample);
//        }
//        assert(dynamic_paras->twiss_ref.bet_x>0&& dynamic_paras->twiss_ref.bet_y>0
//               &&"Need to define the TWISS parameters for the reference point for simulations using the particle model!");
//        ion_beam_model_MonteCarlo_Gaussian(n_sample, ion, dynamic_paras->twiss_ref);
//        if (rdn.get() == nullptr) rdn.reset(new double[n_sample]);
        break;
    }
    default: {
        std::cout<<"Wrong dynamic model!"<<std::endl;
        assert(false);
    }
    }
    return 0;
}

//void ibs_kick(int n_sample, double rate, double twiss, double dt, double emit, double* p) {
//    int k = rate>0?2:-2;
//    double theta = sqrt(k*rate*dt*emit/twiss);
//    gaussian_random(n_sample, rdn.get(), 1, 0);
//    for(int i=0; i<n_sample; ++i) p[i] += theta*rdn[i];
//}

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
        if(ion.bunched()) ion.set_sigma_s(sigma_s);
        //resample the ions
        if(dynamic_paras->ecool()) ion_sample(*ecool_paras, ion, ring, cooler);
        break;
    }
//    case DynamicModel::MODEL_BEAM : { }
    case DynamicModel::PARTICLE : {
//        unsigned int n_sample = n_ion_model;
//        int n_sample = dynamic_paras->n_sample();
//        if(dynamic_paras->ecool()) n_sample = ecool_paras->n_sample();
//        double p0 = ion.gamma()*ion.mass()*1e6*k_e*ion.beta()/k_c;

        if(dynamic_paras->ecool()) {
             double freq = k_c*ion.beta()/ring.circ()*cooler.section_number();
             //restore the coordinates
             restore_cord(t_cooler);
//             for(unsigned int i=0; i<n_sample; ++i) {
//                 xp[i] -= force_x[i]*t_cooler/p0;
//                 yp[i] -= force_y[i]*t_cooler/p0;
//                 dp_p[i] -= force_z[i]*t_cooler/p0;
//             }
             //Adjust the frequency for bunched electron to coasting ion
             if(ebeam.bunched()&&(!ion.bunched())) {
                adjust_freq(freq, ebeam);
//                double sample_length = ebeam.length();
//                double bunch_separate = ecool_paras->bunch_separate();
//                if(sample_length<0) perror("electron bunch length must be positive!");
//                if(bunch_separate>sample_length) {
//                    double duty_factor = sample_length/bunch_separate;
//                    freq *= duty_factor;
//                }
//                else {
//                    perror("Electron bunch length is larger than the distance between electron bunches");
//                }
             }
            //Apply cooling kicks
            if(dynamic_paras->ecool()) {
                apply_cooling_kick(t_cooler, freq, dt);
//                for(unsigned int i=0; i<n_sample; ++i) {
//                    xp[i] = xp[i]!=0?xp[i]*exp(force_x[i]*t_cooler*dt*freq/(xp[i]*p0)):xp[i];
//                    yp[i] = yp[i]!=0?yp[i]*exp(force_y[i]*t_cooler*dt*freq/(yp[i]*p0)):yp[i];
//                    dp_p[i] = dp_p[i]!=0?dp_p[i]*exp(force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0)):dp_p[i];
//                }
            }
         }

        //Apply ibs kicks
        if(dynamic_paras->ibs()) {
            apply_ibs_kick(dt, ion, r_ibs);
//            assert(dynamic_paras->twiss_ref.bet_x>0&& dynamic_paras->twiss_ref.bet_y>0
//                   &&"TWISS parameters for the reference point not defined! Define twiss_ref.");
//            double rx_ibs = r_ibs.at(0);
//            double ry_ibs = r_ibs.at(1);
//            double rs_ibs = r_ibs.at(2);
//
//            ibs_kick(n_sample, rx_ibs, dynamic_paras->twiss_ref.bet_x, dt, ion.emit_x(), xp.get());
//            ibs_kick(n_sample, ry_ibs, dynamic_paras->twiss_ref.bet_y, dt, ion.emit_y(), yp.get());
//            ibs_kick(n_sample, rs_ibs, 1, dt, ion.dp_p()*ion.dp_p(), dp_p.get());

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
//        double emit_x = emit(x_bet.get(), xp_bet.get(), n_sample);
//        double emit_y = emit(y_bet.get(), yp_bet.get(), n_sample);
//        dp = sqrt(emit_p(dp_p.get(), n_sample));
//        if(ion.bunched()) sigma_s = sqrt(emit_p(ds.get(), n_sample));
//
//        ion.set_emit_x(emit_x);
//        ion.set_emit_y(emit_y);
//        ion.set_dp_p(dp);
//        if(ion.bunched()) ion.set_sigma_s(sigma_s);

        break;
    }
    default: {
        assert(false&&"Wrong dynamic model!");
    }
    }
    return 0;
}

void output(double t, std::vector<double> &emit, std::vector<double> &r, bool bunched, std::ofstream &outfile) {
    outfile<<t<<' '<<emit.at(0)<<' '<<emit.at(1)<<' '<<emit.at(2)<<' ';
    if(bunched) outfile<<emit.at(3)<<' ';
    else outfile<<0<<' ';
    outfile<<r.at(0)<<' '<<r.at(1)<<' '<<r.at(2)<<' ';
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
        <<"&column name=rs, type=double, units=1/s, description=\"lonitudinal expansion rate\", &end"<<endl
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

    //Prepare outputting
    std::ofstream outfile;
    outfile.open(dynamic_paras->output_file());
    output_sddshead(n_step+1, outfile);
    outfile.precision(10);
    outfile<<std::showpos;
    outfile<<std::scientific;

    //Set the twiss parameters for the reference point in model beam simulation
    referece_twiss(cooler);

    if(!ecool && dynamic_paras->model()==DynamicModel::PARTICLE)
        assert(dynamic_paras->n_sample()>0 && "Sample particle number must be greater than zero!");

//    if(!ecool && dynamic_paras->model()==DynamicModel::RMS) {}
//    else assert(dynamic_paras->n_sample()>0 && "Sample particle number must be greater than zero!");

//    if(ecool) ecool_paras->set_n_sample(dynamic_paras->n_sample());

    //sample the ion
    sample_the_ions(ion, ring, cooler);

    //Start tracking
    std::cout<<"Start dynamic simulation ... "<<std::endl;
    dynamic_flag = true;

    for(int i=0; i<n_step+1; ++i) {
        if (ion_save_itvl>0 && i%ion_save_itvl==0) save_ions_sdds(dynamic_paras->n_sample(), "ions"+std::to_string(i)+".txt");
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
        if (output_itvl==1) output(t, emit, r, ion.bunched(), outfile);
        else if(i%output_itvl==0) output(t, emit, r, ion.bunched(), outfile);
        if (ion_save_itvl>0 && i%ion_save_itvl==0 && dynamic_paras->model()==DynamicModel::PARTICLE)
            save_ions_sdds(dynamic_paras->n_sample(), "ions"+std::to_string(i)+".txt");

        //Update beam parameters and particles
        update_beam(i, ion, ring, cooler, ebeam, r_ibs, r_ecool);
        t += dt;
        std::cout<<i<<std::endl;
    }

    if (ion_save_itvl>0 && n_step%ion_save_itvl!=0 && dynamic_paras->model()==DynamicModel::PARTICLE)
        save_ions_sdds(dynamic_paras->n_sample(), "ions"+std::to_string(n_step)+".txt");
    dynamic_flag = false;
    outfile.close();
    std::cout<<"Finished dynamic simulation."<<std::endl;

    return 0;
}


