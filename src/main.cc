#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include "dynamic.h"
#include "ecooling.h"
#include "ibs.h"
#include "math_parser.h"
#include "muParserDLL.h"
#include "ring.h"
#include "ui.h"

using std::string;

extern DynamicParas * dynamic_paras;
extern IBSParas * ibs_paras;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern Luminosity *luminosity_paras;
//extern std::unique_ptr<Twiss> twiss_ref;
//extern int n_ion_model;

extern std::map<std::string, Section> sections;
extern muParserHandle_t math_parser;

//extern std::vector<std::string> sections;
//extern std::vector<std::string> ion_pars;


enum class Test {IBS, ECOOL, BOTH, DYNAMICIBS, DYNAMICECOOL, DYNAMICBOTH, DYNAMICIBSBUNCHED};

int main(int argc, char** argv) {

    if(argc>1) {
        std::ifstream input_file(argv[1]);
        string line;
        Section sec_flag = Section::NONE;
        Set_ptrs ptrs;
//        muParserHandle_t hParser = NULL;
//        Set_ion *ion_ptr = nullptr;
        if (math_parser == NULL) {
            initialize_parser(math_parser);
        }
        while (std::getline(input_file,line)) {
            if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
            if (!line.empty()) {
                line = remove_comments(line);
                line = trim_blank(line);
                line = trim_tab(line);
                string line_orgn = line;
                str_toupper(line);
//                    std::cout<<line<<std::endl;
                if (!line.empty()) {
                    if (sections.find(line)!=sections.end()) {
                        sec_flag = sections[line];
                        switch (sec_flag) {
                            case Section::SECTION_ION: {
                                if (ptrs.ion_ptr.get() == nullptr) ptrs.ion_ptr.reset(new Set_ion());
                                break;
                            }
                            case Section::SECTION_RING: {
                                if (ptrs.ring_ptr.get() == nullptr) ptrs.ring_ptr.reset(new Set_ring());
                                break;
                            }
                            case Section::SECTION_IBS: {
                                if (ptrs.ibs_ptr.get() == nullptr) ptrs.ibs_ptr.reset(new Set_ibs());
                                break;
                            }
                            case Section::SECTION_COOLER: {
                                if (ptrs.cooler_ptr.get() == nullptr) ptrs.cooler_ptr.reset(new Set_cooler());
                                break;
                            }
                            case Section::SECTION_E_BEAM: {
                                if (ptrs.e_beam_ptr.get() == nullptr) ptrs.e_beam_ptr.reset(new Set_e_beam());
                                break;
                            }
                            case Section::SECTION_ECOOL: {
                                if (ptrs.ecool_ptr.get() == nullptr) ptrs.ecool_ptr.reset(new Set_ecool());
                                break;
                            }
                            case Section::SECTION_LUMINOSITY: {
                                if (ptrs.luminosity_ptr.get() == nullptr) ptrs.luminosity_ptr.reset(new Set_luminosity());
                                break;
                            }
                            case Section::SECTION_SIMULATION: {
                                if (ptrs.dynamic_ptr.get() == nullptr) ptrs.dynamic_ptr.reset(new Set_dynamic());
                                break;
                            }
                            case Section::SECTION_RUN: {
                                set_section_run(ptrs);
                                break;
                            }
                            case Section::SECTION_SCRATCH: {
//                                if (math_parser == NULL) {
//                                    initialize_parser(math_parser);
//                                }
                                break;
                            }
                            default : {
                                assert(false && "WRONG SECTION NAME!");
                            }
                        }
                    }
                    else {
                        assert(sec_flag>Section::NONE && "Script file must start with a section!");
//                        string::size_type idx = line.find("=");
//                        assert((idx!=string::npos) && "Wrong command!");
//                        string var = line.substr(0, idx);
//                        string val = line.substr(idx+1);
//                        std::cout<<var<<" "<<val<<std::endl;

                        switch (sec_flag) {
                            case Section::SECTION_ION: {
                                define_ion_beam(line, ptrs.ion_ptr.get());
                                break;
                            }
                            case Section::SECTION_RING: {
                                define_ring(line_orgn, ptrs.ring_ptr.get());
                                break;
                            }
                            case Section::SECTION_IBS: {
                                set_ibs(line, ptrs.ibs_ptr.get());
                                break;
                            }
                            case Section::SECTION_COOLER: {
                                define_cooler(line, ptrs.cooler_ptr.get());
                                break;
                            }
                            case Section::SECTION_E_BEAM: {
                                define_e_beam(line_orgn, ptrs.e_beam_ptr.get());
                                break;
                            }
                            case Section::SECTION_ECOOL: {
                                set_ecool(line, ptrs.ecool_ptr.get());
                                break;
                            }
                            case Section::SECTION_LUMINOSITY: {
                                set_luminosity(line, ptrs.luminosity_ptr.get());
                                break;
                            }
                            case Section::SECTION_SIMULATION: {
                                set_simulation(line, ptrs.dynamic_ptr.get());
                                break;
                            }
                            case Section::SECTION_RUN: {
                                run(line, ptrs);
                                break;
                            }
                            case Section::SECTION_SCRATCH: {
                                parse(line, math_parser);
                                break;
                            }
                        }

                    }
                }
            }
        }


    }
    else {
        Test test = Test::BOTH;
        switch (test){
            case Test::BOTH: {
                // define proton beam;
                double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
                int Z;
                Z = 1;
                m0 = 938.272;
                KE = 100e3;
                emit_nx0 = 1.2e-6;
                emit_ny0 = 0.6e-6;
                dp_p0 = 5e-4;
                sigma_s0 = 0.84e-2;
//                sigma_s0 = 2.5e-2;
                N_ptcl = 6.56E9;
                Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

                 // define the lattice of the proton ring
                std::string filename = "MEICColliderRedesign1IP.tfs";
                Lattice lattice(filename);

                //Define the ring
                Ring ring(lattice, p_beam);

                //Set IBS parameters.
                int nu = 100;
                int nv = 100;
                int nz = 40;
                double log_c = 39.9/2;
                ibs_paras = new IBSParas(nu, nv, log_c);

//                //Calculate IBS rate.
//
//                double rx_ibs, ry_ibs, rz_ibs;
//                ibs_rate(lattice, p_beam, *ibs_paras, rx_ibs, ry_ibs, rz_ibs);
//                std::cout<<"ibs rate: "<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;


                //define the cooler
                double cooler_length = 60;
                double n_section = 1;
                double magnetic_field = 1;
                double beta_h = 100;
                double beta_v = 100;
                double dis_h = 0;
//                double dis_h = 2.5;
                double dis_v = 0;
                Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

//////                //define electron beam
//                double n_electron = 2.62E9;
//                double sigma_x = 0.035E-2;
//                double sigma_y = 0.035E-2;
//                double sigma_s = 0.84E-2;
//                GaussianBunch gaussian_bunch(n_electron, sigma_x, sigma_y, sigma_s);
//                double gamma_e = p_beam.gamma();
//                double tmp_tr = 0.1;
//                double tmp_long = 0.5;
//                EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);

                std::string electron_file = "electrons.dat";
                double ns = 1e6;
                double n_electron = 2.62E9;
                int s = 200;
                int line_skip = 0;
                ParticleBunch particle_bunch(n_electron, electron_file);
                particle_bunch.set_s(s);
                particle_bunch.set_skip(line_skip);
                particle_bunch.set_binary(true);
                particle_bunch.load_particle(ns);

//                ParticleBunch particle_bunch(n_electron, electron_file, ns, line_skip, false, 1000, s);
                double gamma_e = p_beam.gamma();
                EBeam e_beam(gamma_e, particle_bunch);


                unsigned int n_sample = 10000;
                ecool_paras = new EcoolRateParas(n_sample);
                //define friction force formula
                force_paras = new ForceParas(ForceFormula::PARKHOMCHUK);

                double rate_x, rate_y, rate_s;
                ecooling_rate(*ecool_paras, *force_paras, p_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
                std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;

//                rate_x += rx_ibs;
//                rate_y += ry_ibs;
//                rate_s += rz_ibs;
//                std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;



    //            double t = 7200;
    //            int n_step = 720;
    //            bool ibs = true;
    //            bool ecool = true;
    //            dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
    //            dynamic_paras->set_model(DynamicModel::MODEL_BEAM);
    //
    //            char file[100] = "tr-0.5eV.txt";
    //            std::ofstream outfile;
    //            outfile.open(file);
    //            dynamic(p_beam, cooler, e_beam, ring, outfile);
    ////            dynamic(p_beam, cooler, e_beam, ring, outfile);
    //            outfile.close();


                break;

            }
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
//                Z = 1;
//                m0 = 938.272;
//                KE = 100e3;
//                emit_nx0 = 6.00332e-006;
//                emit_ny0 = 3.01154e-007;
//                dp_p0 = 0.000997401;
//                sigma_s0 = 0.0284972;
//                N_ptcl = 6.56E9;

                 Z = 1;
                m0 = 938.272;
                KE = 8000;
                emit_nx0 = 2.2e-6;
                emit_ny0 = 2.2e-6;
                dp_p0 = 0.001*0.6;
                N_ptcl = 6.58E11;
                sigma_s0 = 7;
                Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

                 // define the lattice of the proton ring
                std::string filename = "MEICColliderRedesign1IP.tfs";
                Lattice lattice(filename);

                //Define the ring
                Ring ring(lattice, p_beam);

                //Set IBS parameters.
                int nu = 100;
                int nv = 100;
                int nz = 40;
                double log_c = 39.2/2;
                IBSParas ibs_paras(nu, nv, nz);

//                //Calculate IBS rate.
//                std::chrono::steady_clock::time_point start, end;
//                start = std::chrono::steady_clock::now();
//
                double rx_ibs, ry_ibs, rz_ibs;
//                ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);
//
//                end = std::chrono::steady_clock::now();
//                auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//                std::cout<<"IBS 3D integral: "<<t1<<std::endl;
//
//                std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

                ibs_paras.set_log_c(log_c);
                ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);
                std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;


                ibs_paras.set_k(0.2);
                ibs_rate(lattice, p_beam, ibs_paras, rx_ibs, ry_ibs, rz_ibs);
                std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
                break;
            }
            case Test::DYNAMICIBSBUNCHED: {
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
                ibs_paras = new IBSParas(nu, nv, log_c);

                dynamic_paras = new DynamicParas(3600, 360, true, false);

//                char file[100] = "test_dynamic_ibs.txt";
//                std::ofstream outfile;
//                outfile.open(file);
                Cooler *cooler=nullptr;
                EBeam *e_beam=nullptr;
                dynamic(p_beam, *cooler, *e_beam, ring);
//                outfile.close();
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

//                char file[100] = "test_dynamic_ibs.txt";
//                std::ofstream outfile;
//                outfile.open(file);
                Cooler *cooler=nullptr;
                EBeam *e_beam=nullptr;
//                dynamic(p_beam, *cooler, *e_beam, ring, outfile);
//                outfile.close();

                dynamic_paras->set_output_file("test_dynamic_ibs.txt");
                dynamic(p_beam, *cooler, *e_beam, ring);
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

//                char file[100] = "test_dynamic_ecool_DC_model_beam.txt";
//                std::ofstream outfile;
//                outfile.open(file);
//                dynamic(p_beam, cooler, e_beam, ring, outfile);
//                outfile.close();

                dynamic_paras->set_output_file("test_dynamic_ecool_DC_model_beam.txt");
                dynamic(p_beam, cooler, e_beam, ring);

                break;
            }
            case Test::DYNAMICBOTH: {

                srand(time(NULL));
    //            srand(0);

                // define proton beam;
                double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
                int Z;
                Z = 1;
                m0 = 938.272;
                KE = 100e3;
    //            // CM energy 44.7 GeV
    //            emit_nx0 = 0.5e-6;
    //            emit_ny0 = 0.12e-6;
    //            dp_p0 = 0.0008;
    //            N_ptcl = 0.98E10*0.59;
    //            sigma_s0 = 2.5E-2;

                double factor = 2.25;
                // CM energy 44.7 GeV
                emit_nx0 = 0.5e-6*factor;
                emit_ny0 = 0.15e-6*factor;
                dp_p0 = 0.0008;
                N_ptcl = 0.98E10*0.93;
                sigma_s0 = 1.5E-2;

    //            // CM energy 63.3 GeV
    //            emit_nx0 = 1.25e-6;
    //            emit_ny0 = 0.38e-6;
    //            dp_p0 = 0.0008;
    //            N_ptcl = 3.9E10*0.25;
    //            sigma_s0 = 2.5E-2;
                Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);
                std::cout<<"Normalized emittance: "<<p_beam.emit_nx()<<' '<<p_beam.emit_ny()<<std::endl;
                std::cout<<"Geometric emittance: "<<p_beam.emit_x()<<' '<<p_beam.emit_y()<<std::endl;

                // define the lattice of the proton ring
    //            std::string filename = "MEICBoosterRedesign.tfs";
                std::string filename = "MEICColliderRedesign1IP.tfs";
                Lattice lattice(filename);

                //Define the ring
                Ring ring(lattice, p_beam);

    //            //define the cooler
                double cooler_length = 60;
                double n_section = 1;
                double magnetic_field = 1;
                double beta_h = 60;
                double beta_v = 200;
                double dis_h = 2.0;
                double dis_v = 0.6;
                Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);
                std::cout<<"Ion beam size at cooler: "<<sqrt(cooler.beta_h()*p_beam.emit_x())
                        <<' '<<sqrt(cooler.beta_v()*p_beam.emit_y())<<std::endl<<std::endl;

    //            //define electron beam
                double length = 0.02;
                double radius = 0.000528*sqrt(factor);
                std::cout<<"Electron beam radius: "<<radius<<std::endl;
                double q_e = 2.0e-9;
                double current = q_e*p_beam.beta()*k_c/length;
                UniformBunch uniform_bunch(current, radius, length);
                double gamma_e = p_beam.gamma();
                double tmp_tr = 0.246;
                double tmp_long = 0.184;
                EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_bunch);

    //            double n_e = 1.248E10;
    //            double sigma_e_x = 0.035E-2;
    //            double sigma_e_y = 0.035E-2;
    //            double sigma_e_s = 0.84E-2;
    //            GaussianBunch gaussian_bunch(n_e, sigma_e_x, sigma_e_y, sigma_e_s);
    //            double gamma_e = p_beam.gamma();
    //            double tmp_tr = 0.5;
    //            double tmp_long = 0.1;
    //            EBeam e_beam(gamma_e, tmp_tr, tmp_long, gaussian_bunch);
    //
                //define cooling model: monte carlo
                unsigned int n_sample = 100000;
                ecool_paras = new EcoolRateParas(n_sample);
    //            //define friction force formula
                force_paras = new ForceParas(ForceFormula::PARKHOMCHUK);

                double rate_x, rate_y, rate_s;
                ecooling_rate(*ecool_paras, *force_paras, p_beam, cooler, e_beam, ring, rate_x, rate_y, rate_s);
                std::cout<<"cooling rate: [1/s] "<<rate_x<<" "<<rate_y<<" "<<rate_s<<std::endl;

                //Set IBS parameters.
                int nu = 100;
                int nv = 100;
                int nz = 40;
                double log_c = 40.4/2;      //100 GeV, 63.3 GeV CM Energy
    //            double log_c = 39.2/2;    //100 GeV
                ibs_paras = new IBSParas(nu, nv, log_c);
    //            ibs_paras = new IBSParas(nu, nv, nz);
                ibs_paras->set_k(0.4);

                double rx_ibs, ry_ibs, rz_ibs;
                ibs_rate(lattice, p_beam, *ibs_paras, rx_ibs, ry_ibs, rz_ibs);
                std::cout<<"IBS rate: [1/s] "<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
                std::cout<<"Total rate: [1/s] "<<rx_ibs+rate_x<<' '<<ry_ibs+rate_y<<' '<<rz_ibs+rate_s<<std::endl<<std::endl;

//                return 0;

                //define dynamic simulation
                double t = 1200*2;
                int n_step = 600*2;
                bool ibs = true;
                bool ecool = true;
                dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
                dynamic_paras->set_model(DynamicModel::MODEL_BEAM);
//                dynamic_paras->set_model(DynamicModel::RMS);

//                char file[100] = "Collider_100GeV_strong_cooling_baseline_44Ecm_2nC_rms_02_long_cool_compensated.txt";
//                std::ofstream outfile;
//                outfile.open(file);
//                Cooler *cooler;
//                EBeam *e_beam;
//                if(twiss_ref.get()==nullptr) twiss_ref.reset(new Twiss());
//                twiss_ref->bet_x = 10;
//                twiss_ref->bet_y = 10;
//                twiss_ref->disp_x = 5;
                dynamic_paras->twiss_ref.bet_x = 10;
                dynamic_paras->twiss_ref.bet_y = 10;
                dynamic_paras->twiss_ref.disp_x = 5;
//                dynamic_paras->set_n_sample(5000);
//                dynamic(p_beam, *cooler, *e_beam, ring, outfile);
//                dynamic(p_beam, cooler, e_beam, ring, outfile);
//                outfile.close();

                dynamic_paras->set_output_file("Collider_100GeV_strong_cooling_baseline_44Ecm_2nC_rms_02_long_cool_compensated.txt");
                dynamic(p_beam, cooler, e_beam, ring);

                break;
            }
            default: {
                assert(false);
            }
        }
    }

    //Pause the system
    std::cout<<std::endl<<"Press the Enter key to close the window."<<std::endl;
    std::cin.ignore( std::numeric_limits< std::streamsize >::max( ), '\n' );

    return 0;
}

