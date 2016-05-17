#include <chrono>
#include<fstream>
#include "dynamic.h"
#include "ecooling.h"
#include "ibs.h"
#include "ring.h"

extern DynamicParas * dynamic_paras;
extern IBSParas * ibs_paras;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern Twiss * twiss_ref;
extern int n_ion_model;

enum class Test {IBS, ECOOL, BOTH, DYNAMICIBS, DYNAMICECOOL, DYNAMICBOTH};

int main() {

    Test test = Test::DYNAMICBOTH;
    switch (test){
        case Test::ECOOL: {
            //********************************
            // Test Electron Cooling Rate
            //********************************

            //define coasting 12C6+ beam
            int n_charge = 6;
            double n_mass = 12;
            double kinetic_energy = 30*n_mass;
            double gamma = 1+kinetic_energy/(n_mass*k_u);
            double beta = sqrt(1-1/(gamma*gamma));
            double emit_nx0 = beta*gamma*5e-6;
            double emit_ny0 = beta*gamma*5e-6;
            double dp_p0 = 0.0004;
            double n_ptcl = 5E8;
            Beam c_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, n_ptcl);

            //define the ring
            std::string filename = "csrm.tfs";
            Lattice lattice(filename);
            Ring ring(lattice, c_beam);

            //define the cooler
            double cooler_length = 3.4;
            double n_section = 1;
            double magnetic_field = 0.039;
            double beta_h = 10;
            double beta_v = 17;
            double dis_h = 0;
            double dis_v = 0;
            Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

            //define electron beam
            double current = 0.03;
            double radius = 0.025;
            double neutralisation = 0;
            UniformCylinder uniform_cylinder(current, radius, neutralisation);
            double gamma_e = c_beam.gamma();
            double tmp_tr = 0.05;
            double tmp_long = 0.1;
            EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_cylinder);


            //define cooling model
//            unsigned int n_sample = 5000;
//            EcoolRateParas ecool_rate_paras(n_sample);
            unsigned int n_tr = 100;
            unsigned int n_long = 100;
            EcoolRateParas ecool_rate_paras(n_tr, n_long);

            ForceParas force_paras(ForceFormula::PARKHOMCHUK);
            double rate_x, rate_y, rate_s;
            ecooling_rate(ecool_rate_paras, force_paras, c_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
            std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;

            break;
        }
        case Test::IBS: {
            //********************************
            // Test IBS rate
            //********************************
            // define proton beam;
            double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
            int Z;
            Z = 1;
            m0 = 938.272;
            KE = 3e4;
            emit_nx0 = 0.4962695094e-6;
            emit_ny0 = 0.4962695094e-6;
            dp_p0 = 4e-4;
            sigma_s0 = 1.994525702e-2;
            N_ptcl = 6.56E9;
            Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

             // define the lattice of the proton ring
            std::string filename = "MEICColliderRedesign1IP.tfs";
            Lattice lattice(filename);

            //Define the ring
            Ring ring(lattice, p_beam);

            //Set IBS parameters.
            int nu = 200;
            int nv = 200;
            int nz = 40;
            double log_c = 39.9/2;
            IBSParas ibs_paras(nu, nv, nz);

            //Calculate IBS rate.
            std::chrono::steady_clock::time_point start, end;
            start = std::chrono::steady_clock::now();

            double rx_ibs, ry_ibs, rz_ibs;
            config_ibs(lattice);
            ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);

            end = std::chrono::steady_clock::now();
            auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            std::cout<<"IBS 3D integral: "<<t1<<std::endl;

            std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;


            ibs_paras.set_k(1);
            ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);
            std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

            ibs_paras.set_log_c(log_c);
            ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);
            std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
            end_ibs();
            break;
        }
        case Test::DYNAMICIBS : {
            // define proton beam;
            double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
            int Z;
            Z = 1;
            m0 = 938.272;
            KE = 800;
            emit_nx0 = 1.039757508e-6;
            emit_ny0 = 1.039757508e-6;
            dp_p0 = 2e-3;
            N_ptcl = 3.6E11;
            Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, N_ptcl);

             // define the lattice of the proton ring
            std::string filename = "MEICBoosterRedesign.tfs";
            Lattice lattice(filename);

            //Define the ring
            Ring ring(lattice, p_beam);

            //Set IBS parameters.
            int nu = 200;
            int nv = 200;
            int nz = 40;
            double log_c = 44.8/2;
            ibs_paras = new IBSParas(nu, nv, log_c);
            ibs_paras->set_k(1.0);

            dynamic_paras = new DynamicParas(3600, 360, true, false);

            char file[100] = "test_dynamic_ibs.txt";
            std::ofstream outfile;
            outfile.open(file);
            Cooler *cooler=nullptr;
            EBeam *e_beam=nullptr;
            dynamic(p_beam, *cooler, *e_beam, ring, outfile);
            outfile.close();
            break;
        }
        case Test::DYNAMICECOOL: {
            // define proton beam;
            double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
            int Z;
            Z = 1;
            m0 = 938.272;
            KE = 800;
            emit_nx0 = 1.039757508e-6;
            emit_ny0 = 1.039757508e-6;
            dp_p0 = 0.002012615391;
            N_ptcl = 3.6E11;
            Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, N_ptcl);

            // define the lattice of the proton ring
            std::string filename = "MEICBoosterRedesign.tfs";
            Lattice lattice(filename);

            //Define the ring
            Ring ring(lattice, p_beam);

            //define the cooler
            double cooler_length = 10;
            double n_section = 1;
            double magnetic_field = 0.1;
            double beta_h = 10;
            double beta_v = 10;
            double dis_h = 0;
            double dis_v = 0;
            Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

            //define electron beam
            double current = 2;
            double radius = 0.008;
            double neutralisation = 0;
            UniformCylinder uniform_cylinder(current, radius, neutralisation);
            double gamma_e = p_beam.gamma();
            double tmp_tr = 0.1;
            double tmp_long = 0.1;
            EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_cylinder);

//             //define cooling model: single particle
//            unsigned int n_tr = 100;
//            unsigned int n_long = 100;
//            ecool_paras = new EcoolRateParas(n_tr, n_long);
            //define cooling model: monte carlo
            unsigned int n_sample = 40000;
            ecool_paras = new EcoolRateParas(n_sample);
            //define friction force formula
            force_paras = new ForceParas(ForceFormula::PARKHOMCHUK);
            //define dynamic simulation
            dynamic_paras = new DynamicParas(60, 120, false, true);
            dynamic_paras->set_model(DynamicModel::MODEL_BEAM);

            char file[100] = "test_dynamic_ecool_DC_model_beam.txt";
            std::ofstream outfile;
            outfile.open(file);
            dynamic(p_beam, cooler, e_beam, ring, outfile);
            outfile.close();

            break;
        }
        case Test::DYNAMICBOTH: {
            // define proton beam;
            double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
            int Z;
            Z = 1;
            m0 = 938.272;
            KE = 30000;
            emit_nx0 = 0.4515633419e-6;
            emit_ny0 = 0.4515633419e-6;
            dp_p0 = 0.0004957205723;
            N_ptcl = 6.56E9;
            sigma_s0 = 2.509282121E-2;
            Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

            // define the lattice of the proton ring
//            std::string filename = "MEICBoosterRedesign.tfs";
            std::string filename = "MEICColliderRedesign1IP.tfs";
            Lattice lattice(filename);

            //Define the ring
            Ring ring(lattice, p_beam);

//            //define the cooler
//            double cooler_length = 30;
//            double n_section = 2;
//            double magnetic_field = 1;
//            double beta_h = 100;
//            double beta_v = 100;
//            double dis_h = 5;
//            double dis_v = 0;
//            Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);
//
//            //define electron beam
//            double n_e = 1.25E10;
//            double sigma_e_x = 0.0012;
//            double sigma_e_y = 0.0012;
//            double sigma_e_s = 0.025;
//            GaussianBunch gaussian_bunch(n_e, sigma_e_x, sigma_e_y, sigma_e_s);
//            double gamma_e = p_beam.gamma();
//            double tmp_tr = 0.5;
//            double tmp_long = 0.1;
//            EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);
//
//             //define cooling model: single particle
////            unsigned int n_tr = 100;
////            unsigned int n_long = 100;
////            ecool_paras = new EcoolRateParas(n_tr, n_long);
//            //define cooling model: monte carlo
//            unsigned int n_sample = 10000;
//            ecool_paras = new EcoolRateParas(n_sample);
//            //define friction force formula
//            force_paras = new ForceParas(ForceFormula::PARKHOMCHUK);
            //define dynamic simulation
            double t = 3600;
            int n_step = 360;
            bool ibs = true;
            bool ecool = false;
            dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
            dynamic_paras->set_model(DynamicModel::MODEL_BEAM);

            //Set IBS parameters.
            int nu = 100;
            int nv = 100;
            int nz = 40;
            double log_c = 39.2/2;
            ibs_paras = new IBSParas(nu, nv, log_c);

            char file[100] = "test_dynamic_ibs_bunched_model_beam.txt";
            std::ofstream outfile;
            outfile.open(file);
            Cooler *cooler;
            EBeam *e_beam;
            twiss_ref = new Twiss();
            twiss_ref->bet_x_ = 10;
            twiss_ref->bet_y_ = 10;
//            twiss_ref->disp_x_ = 5;
            n_ion_model = 5000;
            dynamic(p_beam, *cooler, *e_beam, ring, outfile);
//            dynamic(p_beam, cooler, e_beam, ring, outfile);
            outfile.close();

            break;
        }
        default: {
            assert(false);
        }
    }

    return 0;
}

