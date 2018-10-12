#include "ui.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "beam.h"
#include "constants.h"
#include "ibs.h"

using std::string;

extern DynamicParas *dynamic_paras;
extern IBSParas *ibs_paras;
extern EcoolRateParas *ecool_paras;
extern ForceParas *force_paras;

muParserHandle_t math_parser = NULL;
//std::vector<std::string> sections = {"SECTION_ION", "SECTION_RING", "SECTION_COOLER"};
std::vector<string> ION_ARGS = {"CHARGE_NUMBER", "MASS", "KINETIC_ENERGY", "NORM_EMIT_X", "NORM_EMIT_Y",
    "MOMENTUM_SPREAD", "PARTICLE_NUMBER", "RMS_BUNCH_LENGTH"};
std::vector<string> RUN_COMMANDS = {"CREATE_ION_BEAM", "CREATE_RING", "CREATE_E_BEAM", "CREATE_COOLER",
    "CALCULATE_IBS", "CALCULATE_ECOOL", "TOTAL_EXPANSION_RATE", "RUN_SIMULATION"};
std::vector<string> RING_ARGS = {"LATTICE", "QX", "QY", "QS", "GAMMA_TR", "RF_V", "RF_H", "RF_PHI"};
std::vector<string> IBS_ARGS = {"NU","NV","NZ","LOG_C","COUPLING"};
std::vector<string> COOLER_ARGS = {"LENGTH", "SECTION_NUMBER", "MAGNETIC_FIELD", "BET_X", "BET_Y", "DISP_X", "DISP_Y",
    "ALPHA_X", "ALPHA_Y", "DISP_DX", "DISP_DY"};
//std::vector<string> SCRATCH_COMMANDS = {"PRINT", "LIST_VAR", "LIST_CONST"};
//std::vector<string> E_BEAM_SHAPE_ARGS = {"SHAPE", "RADIUS", "CURRENT", "SIMGA_X", "SIGMA_Y", "SIGMA_Z", "LENGTH", "E_NUMBER"};
std::vector<string> E_BEAM_SHAPE_TYPES = {"DC_UNIFORM", "BUNCHED_GAUSSIAN", "BUNCHED_UNIFORM", "BUNCHED_UNIFORM_ELLIPTIC",
    "DC_UNIFORM_HOLLOW", "BUNCHED_UNIFORM_HOLLOW", "BUNCHED_USER_DEFINED"};
//std::vector<string> E_BEAM_ARGS = {"GAMMA", "TMP_TR", "TMP_L"};
std::vector<string> E_BEAM_ARGS = {"GAMMA", "TMP_TR", "TMP_L", "SHAPE", "RADIUS", "CURRENT", "SIGMA_X", "SIGMA_Y",
    "SIGMA_Z", "LENGTH", "E_NUMBER", "RH", "RV", "R_INNER", "R_OUTTER", "PARTICLE_FILE", "TOTAL_PARTICLE_NUMBER",
    "BOX_PARTICLE_NUMBER", "LINE_SKIP", "VEL_POS_CORR"};
std::vector<string> ECOOL_ARGS = {"SAMPLE_NUMBER", "FORCE_FORMULA"};
std::vector<string> FRICTION_FORCE_FORMULA = {"PARKHOMCHUK"};
std::vector<string> SIMULATION_ARGS = {"TIME", "STEP_NUMBER", "SAMPLE_NUMBER", "IBS", "E_COOL", "OUTPUT_INTERVAL",
    "SAVE_PARTICLE_INTERVAL", "OUTPUT_FILE", "MODEL", "REF_BET_X", "REF_BET_Y", "REF_ALF_X", "REF_ALF_Y",
    "REF_DISP_X", "REF_DISP_Y", "REF_DISP_DX", "REF_DISP_DY"};
std::vector<string> DYNAMIC_VALUE = {"RMS", "PARTICLE", "MODEL_BEAM", "TURN_BY_TURN"};

std::map<std::string, Section> sections{
    {"SECTION_ION",Section::SECTION_ION},
    {"SECTION_RING",Section::SECTION_RING},
    {"SECTION_COOLER",Section::SECTION_COOLER},
    {"SECTION_RUN",Section::SECTION_RUN},
    {"SECTION_IBS",Section::SECTION_IBS},
    {"SECTION_SCRATCH", Section::SECTION_SCRATCH},
    {"SECTION_E_BEAM", Section::SECTION_E_BEAM},
    {"SECTION_ECOOL", Section::SECTION_ECOOL},
    {"SECTION_SIMULATION",Section::SECTION_SIMULATION}
};

// Remove everything from the first "#" in the string
std::string remove_comments(std::string input_line) {
    return input_line.substr(0,input_line.find('#'));
}

// Trim the spaces at the head and the tail of the string
std::string trim_blank(std::string input_line) {
    using std::string;
    if (input_line.empty()) return input_line;
    string::size_type st = input_line.find_first_not_of(" ");
    if (st == string::npos) return "";
    string::size_type fi = input_line.find_last_not_of(" ")+1;
    return input_line.substr(st, fi);
}


std::string trim_tab(std::string input_line) {
    using std::string;
    if (input_line.empty()) return input_line;
    string::size_type st = input_line.find_first_not_of("\t");
    if (st == string::npos) return "";
    string::size_type fi = input_line.find_last_not_of("\t")+1;
    return input_line.substr(st, fi);
}

void str_toupper(std::string &str) {
    for (auto & c: str) c = toupper(c);
}

void define_e_beam(string &str, Set_e_beam *e_beam_args) {
    assert(e_beam_args!=nullptr && "SECTION_E_BEAM MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_E_BEAM!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(E_BEAM_ARGS.begin(),E_BEAM_ARGS.end(),var)!=E_BEAM_ARGS.end()
           && "WRONG COMMANDS IN SECTION_E_BEAM!");
    if (var == "SHAPE") {
        assert(std::find(E_BEAM_SHAPE_TYPES.begin(),E_BEAM_SHAPE_TYPES.end(),val)!=E_BEAM_SHAPE_TYPES.end()
           && "UNDEFINED ELECTRON BEAM SHAPE!");
        e_beam_args->shape = val;
    }
    else if (var == "PARTICLE_FILE") {
        e_beam_args->particle_file = val;
    }
    else {
        if (math_parser == NULL) {
            if(var == "E_NUMBER") {
                e_beam_args->n = std::stod(val);
            }
            else if (var == "RADIUS") {
                e_beam_args->radius = std::stod(val);
            }
            else if (var == "CURRENT") {
                e_beam_args->current = std::stod(val);
            }
            else if (var == "SIMGA_X") {
                e_beam_args->sigma_x = std::stod(val);
            }
            else if (var == "SIGMA_Y") {
                e_beam_args->sigma_y = std::stod(val);
            }
            else if (var == "SIGMA_Z") {
                e_beam_args->sigma_z = std::stod(val);
            }
            else if (var == "RH") {
                e_beam_args->rh = std::stod(val);
            }
            else if (var == "RV") {
                e_beam_args->rv = std::stod(val);
            }
            else if (var == "LENGTH") {
                e_beam_args->length = std::stod(val);
            }
            else if (var == "GAMMA") {
                e_beam_args->gamma = std::stod(val);
            }
            else if (var == "TMP_TR") {
                e_beam_args->tmp_tr = std::stod(val);
            }
            else if (var == "TMP_L") {
                e_beam_args->tmp_l = std::stod(val);
            }
            else if (var == "R_INNER") {
                e_beam_args->r_inner = std::stod(val);
            }
            else if (var == "R_OUTTER") {
                e_beam_args->r_outter = std::stod(val);
            }
            else if (var == "TOTAL_PARTICLE_NUMBER") {
                e_beam_args->n_particle = std::stoi(val);
            }
            else if (var == "LINE_SKIP") {
                e_beam_args->line_skip = std::stoi(val);
            }
            else if (var == "BOX_PARTICLE_NUMBER") {
                e_beam_args->particle_perbox = std::stoi(val);
            }
            else if (var == "BUFFER_SIZE") {
                e_beam_args->buffer = std::stoi(val);
            }
            else if (var == "VEL_POS_CORR") {
                int v = std::stoi(val);
                switch (v) {
                case 0 : {
                    e_beam_args->corr = false;
                }
                case 1 : {
                    e_beam_args->corr = true;
                }
                default : {
                    assert(false&& "WRONG VALUE FOR VEL_POS_CORR FOR E_BEAM!");
                }
                }
            }
            else if (var == "BINARY_FILE") {
                int v = std::stoi(val);
                switch (v) {
                case 0 : {
                    e_beam_args->binary = false;
                }
                case 1 : {
                    e_beam_args->binary = true;
                }
                default : {
                    assert(false&& "WRONG VALUE FOR BINARY_FILE FOR E_BEAM!");
                }
                }
            }
            else {
                assert(false&&"Wrong arguments in section_e_beam!");
            }
        }
        else {
            mupSetExpr(math_parser, val.c_str());
            if(var == "E_NUMBER") {
                e_beam_args->n = mupEval(math_parser);
            }
            else if (var == "RADIUS") {
                e_beam_args->radius = mupEval(math_parser);
            }
            else if (var == "CURRENT") {
                e_beam_args->current = mupEval(math_parser);
            }
            else if (var == "SIGMA_X") {
                e_beam_args->sigma_x = mupEval(math_parser);
            }
            else if (var == "SIGMA_Y") {
                e_beam_args->sigma_y = mupEval(math_parser);
            }
            else if (var == "SIGMA_Z") {
                e_beam_args->sigma_z = mupEval(math_parser);
            }
            else if (var == "RH") {
                e_beam_args->rh = mupEval(math_parser);
            }
            else if (var == "RV") {
                e_beam_args->rv = mupEval(math_parser);
            }
            else if (var == "LENGTH") {
                e_beam_args->length = mupEval(math_parser);
            }
            else if (var == "GAMMA") {
                e_beam_args->gamma = mupEval(math_parser);
            }
            else if (var == "TMP_TR") {
                e_beam_args->tmp_tr = mupEval(math_parser);
            }
            else if (var == "TMP_L") {
                e_beam_args->tmp_l = mupEval(math_parser);
            }
            else if (var == "R_INNER") {
                e_beam_args->r_inner = mupEval(math_parser);
            }
            else if (var == "R_OUTTER") {
                e_beam_args->r_outter = mupEval(math_parser);
            }
            else if (var == "TOTAL_PARTICLE_NUMBER") {
                e_beam_args->n_particle = mupEval(math_parser);
            }
            else if (var == "LINE_SKIP") {
                e_beam_args->line_skip = mupEval(math_parser);
            }
            else if (var == "BOX_PARTICLE_NUMBER") {
                e_beam_args->particle_perbox = mupEval(math_parser);
            }
            else if (var == "BUFFER_SIZE") {
                e_beam_args->buffer = mupEval(math_parser);
            }
            else if (var == "VEL_POS_CORR") {
                int v = mupEval(math_parser);
                switch (v) {
                case 0 : {
                    e_beam_args->corr = false;
                    break;
                }
                case 1 : {
                    e_beam_args->corr = true;
                    break;
                }
                default : {
                    assert(false&& "WRONG VALUE FOR VEL_POS_CORR FOR E_BEAM!");
                    break;
                }
                }
            }
            else if (var == "BINARY_FILE") {
                int v = mupEval(math_parser);
                switch (v) {
                case 0 : {
                    e_beam_args->binary = false;
                    break;
                }
                case 1 : {
                    e_beam_args->binary = true;
                    break;
                }
                default : {
                    assert(false&& "WRONG VALUE FOR BINARY_FILE FOR E_BEAM!");
                    break;
                }
                }
            }
            else {
                assert(false&&"Wrong arguments in section_e_beam!");
            }
        }
    }
}

void create_e_beam(Set_ptrs &ptrs) {
    assert(ptrs.e_beam_ptr.get()!=nullptr && "MUST DEFINE THE ELECTRON BEAM BEFORE CREATE THE ELECTRON BEAM!");
    std::string shape = ptrs.e_beam_ptr->shape;
    assert(std::find(E_BEAM_SHAPE_TYPES.begin(),E_BEAM_SHAPE_TYPES.end(),shape)!=E_BEAM_SHAPE_TYPES.end()
           && "WRONG ELECTRON BEAM SHAPE!");
    double gamma = ptrs.e_beam_ptr->gamma;
    double tmp_tr = ptrs.e_beam_ptr->tmp_tr;
    double tmp_l = ptrs.e_beam_ptr->tmp_l;
    assert(gamma>0 && tmp_tr >= 0 && tmp_l >= 0 && "WRONG PARAMETER VALUE FOR ELECTRON BEAM!");

    if (shape == "DC_UNIFORM") {
        double current = ptrs.e_beam_ptr->current;
        double radius = ptrs.e_beam_ptr->radius;
        assert(current >= 0 && radius > 0 && "WRONG PARAMETER VALUE FOR DC_UNIFORM SHAPE");
        ptrs.e_beam_shape.reset(new UniformCylinder(current, radius));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "BUNCHED_GAUSSIAN") {
        double n = ptrs.e_beam_ptr->n;
        double sigma_x = ptrs.e_beam_ptr->sigma_x;
        double sigma_y = ptrs.e_beam_ptr->sigma_y;
        double sigma_z = ptrs.e_beam_ptr->sigma_z;
        assert(sigma_x > 0 && sigma_y > 0 && sigma_z > 0 && n > 0 && "WRONG PARAMETER VALUE FOR BUNCHED_GAUSSIAN SHAPE");
        ptrs.e_beam_shape.reset(new GaussianBunch(n, sigma_x, sigma_y, sigma_z));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "BUNCHED_UNIFORM") {
        double current = ptrs.e_beam_ptr->current;
        double radius = ptrs.e_beam_ptr->radius;
        double length = ptrs.e_beam_ptr->length;
        assert(current >= 0 && radius > 0 && length > 0 && "WRONG PARAMETER VALUE FOR BUNCHED_UNIFORM SHAPE");
        ptrs.e_beam_shape.reset(new UniformBunch(current, radius, length));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "BUNCHED_UNIFORM_ELLIPTIC") {
        double current = ptrs.e_beam_ptr->current;
        double rh = ptrs.e_beam_ptr->rh;
        double rv = ptrs.e_beam_ptr->rv;
        double length = ptrs.e_beam_ptr->length;
        assert(current >= 0 && rh > 0 && rv > 0 && length > 0 && "WRONG PARAMETER VALUE FOR BUNCHED_UNIFORM_ELLIPTIC SHAPE");
        ptrs.e_beam_shape.reset(new EllipticUniformBunch(current, rh, rv, length));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "DC_UNIFORM_HOLLOW") {
        double current = ptrs.e_beam_ptr->current;
        double r_inner = ptrs.e_beam_ptr->r_inner;
        double r_outter = ptrs.e_beam_ptr->r_outter;
        assert(r_inner>0 && r_outter>0 && current>=0 && r_outter>r_inner && "WRONG PARAMETER VALUE FOR DC_UNIFORM_HOLLOW SHAPE");
        ptrs.e_beam_shape.reset(new UniformHollow(current, r_inner, r_outter));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "BUNCHED_UNIFORM_HOLLOW") {
        double current = ptrs.e_beam_ptr->current;
        double r_inner = ptrs.e_beam_ptr->r_inner;
        double r_outter = ptrs.e_beam_ptr->r_outter;
        double length = ptrs.e_beam_ptr->length;
        assert(r_inner>0 && r_outter>0 && current>=0 && r_outter>r_inner && length>0 && "WRONG PARAMETER VALUE FOR DC_UNIFORM_HOLLOW SHAPE");
        ptrs.e_beam_shape.reset(new UniformHollowBunch(current, r_inner, r_outter, length));
        ptrs.e_beam.reset(new EBeam(gamma, tmp_tr, tmp_l, *ptrs.e_beam_shape.get()));
    }
    else if(shape == "BUNCHED_USER_DEFINED") {
        double n_electron = ptrs.e_beam_ptr->n;
        std::string filename = ptrs.e_beam_ptr->particle_file;
        int line_skip = ptrs.e_beam_ptr->line_skip;
        int n_particle = ptrs.e_beam_ptr->n_particle;
        int s = ptrs.e_beam_ptr->particle_perbox;
        int buffer = ptrs.e_beam_ptr->buffer;
        double length = ptrs.e_beam_ptr->length;
        assert(n_electron>0 && line_skip>=0 && n_particle>=0 && s>0 && length>=0 && buffer>0 && "WRONG PARAMETER VALUE FOR BUNCHED_USER_DEFINED SHAPE");
        if(length>0)
            ptrs.e_beam_shape.reset(new ParticleBunch(n_electron, filename, n_particle, length, line_skip, s));
        else
            ptrs.e_beam_shape.reset(new ParticleBunch(n_electron, filename, n_particle, line_skip, s));
        ParticleBunch* prtl_bunch = nullptr;
        prtl_bunch = dynamic_cast<ParticleBunch*>(ptrs.e_beam_shape.get());
        if(ptrs.e_beam_ptr->binary)
            prtl_bunch->set_binary(ptrs.e_beam_ptr->binary);
        prtl_bunch->set_s(s);
        prtl_bunch->set_skip(line_skip);
        if(n_particle>0)
            prtl_bunch->load_particle(n_particle);
        else
            prtl_bunch->load_particle();
        if(ptrs.e_beam_ptr->corr) {
            prtl_bunch->set_corr(true);
        }
        ptrs.e_beam.reset(new EBeam(gamma,*ptrs.e_beam_shape.get()));
    }
    std::cout<<"Electron beam created!"<<std::endl;
}

void define_ion_beam(std::string &str, Set_ion *ion_args){
    assert(ion_args!=nullptr && "SECTION_ION MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_ION!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(ION_ARGS.begin(),ION_ARGS.end(),var)!=ION_ARGS.end() && "WRONG COMMANDS IN SECTION_ION!");
//    std::cout<<var<<" "<<val<<std::endl;
    if(math_parser==NULL) {
        if (var=="CHARGE_NUMBER") {
            ion_args->n_charge = std::stoi(val);
        }
        else if (var=="MASS") {
            ion_args->mass = std::stod(val);
        }
        else if (var=="KINETIC_ENERGY") {
            ion_args->k_energy = std::stod(val);
        }
        else if (var=="NORM_EMIT_X") {
            ion_args->emit_nx = std::stod(val);
        }
        else if (var=="NORM_EMIT_Y") {
            ion_args->emit_ny = std::stod(val);
        }
        else if (var=="MOMENTUM_SPREAD") {
            ion_args->dp_p = std::stod(val);
        }
        else if (var=="PARTICLE_NUMBER") {
            ion_args->n_ptcl = std::stod(val);
        }
        else if (var=="RMS_BUNCH_LENGTH") {
            ion_args->ds = std::stod(val);
        }
        else {
            assert(false&&"Wrong arguments in section_ion!");
        }
    }
    else {
        mupSetExpr(math_parser, val.c_str());
        if (var=="CHARGE_NUMBER") {
            ion_args->n_charge = static_cast<int>(mupEval(math_parser));
        }
        else if (var=="MASS") {
            ion_args->mass = mupEval(math_parser);
        }
        else if (var=="KINETIC_ENERGY") {
            ion_args->k_energy = mupEval(math_parser);
        }
        else if (var=="NORM_EMIT_X") {
            ion_args->emit_nx = mupEval(math_parser);
        }
        else if (var=="NORM_EMIT_Y") {
            ion_args->emit_ny = mupEval(math_parser);
        }
        else if (var=="MOMENTUM_SPREAD") {
            ion_args->dp_p = mupEval(math_parser);
        }
        else if (var=="PARTICLE_NUMBER") {
            ion_args->n_ptcl = mupEval(math_parser);
        }
        else if (var=="RMS_BUNCH_LENGTH") {
            ion_args->ds = mupEval(math_parser);
        }
        else {
            assert(false&&"Wrong arguments in section_ion!");
        }
    }
}

void create_ion_beam(Set_ptrs &ptrs) {
    int n_charge = ptrs.ion_ptr->n_charge;
    double mass = ptrs.ion_ptr->mass;
    double k_energy = ptrs.ion_ptr->k_energy;
    double emit_nx = ptrs.ion_ptr->emit_nx;
    double emit_ny = ptrs.ion_ptr->emit_ny;
    double dp_p = ptrs.ion_ptr->dp_p;
    double n_ptcl = ptrs.ion_ptr->n_ptcl;
    double ds = ptrs.ion_ptr->ds;
    assert(n_charge>0 && mass>0 && k_energy>0 && emit_nx>0 && emit_ny>0 && dp_p>=0 && n_ptcl>0
           && "WRONG PARAMETER VALUE FOR ION BEAM!");
    if (ds>0) {
        ptrs.ion_beam.reset(new Beam(n_charge, mass/k_u, k_energy, emit_nx, emit_ny, dp_p, ds, n_ptcl));
        std::cout<<"Bunched ion beam created!"<<std::endl;
    }
    else {
        ptrs.ion_beam.reset(new Beam(n_charge, mass/k_u, k_energy, emit_nx, emit_ny, dp_p, n_ptcl));
        std::cout<<"Coasting ion beam created!"<<std::endl;
    }
}

void create_ring(Set_ptrs &ptrs) {
    string lattice_file = ptrs.ring_ptr->lattice_file;
    assert(!lattice_file.empty() && "WRONG PARAMETER VALUE FOR RING");
    ptrs.lattice.reset(new Lattice(lattice_file));
//    ptrs.lattice = std::make_shared<Lattice>(lattice_file);
//    std::cout<< ptrs.lattice->betx(0) <<std::endl;
    assert(ptrs.ion_beam.get()!=nullptr && "MUST DEFINE THE ION BEFORE CREATE THE RING!");
    ptrs.ring.reset(new Ring(*ptrs.lattice, *ptrs.ion_beam));

    ptrs.tunes.reset(new Tunes());
    ptrs.ring->tunes = ptrs.tunes.get();
    ptrs.rf.reset(new RF());
    ptrs.ring->rf = ptrs.rf.get();
    if(ptrs.ring_ptr->qx>0 || ptrs.ring_ptr->qy>0 || ptrs.ring_ptr->qs>0) {
        assert(ptrs.ring_ptr->qx>0 && ptrs.ring_ptr->qy>0 && "Transverse tunes must be greater than zero!");
//        ptrs.ring.qx = ptrs.ring_ptr->qx;
//        ptrs.ring.qy = ptrs.ring_ptr->qy;
//        ptrs.ring.qs = ptrs.ring_ptr->qs;
//        ptrs.tunes.reset(new Tunes());
        ptrs.tunes->qx = ptrs.ring_ptr->qx;
        ptrs.tunes->qy = ptrs.ring_ptr->qy;
        ptrs.tunes->qs = ptrs.ring_ptr->qs;
//        ptrs.ring->tunes = ptrs.tunes.get();
    }
    if(ptrs.ring_ptr->rf_v>0) {
        assert(ptrs.ring_ptr->rf_v>0 && ptrs.ring_ptr->gamma_tr>0 && "RF Voltage and transition gamma should be greater than zero");
//        ptrs.ring.v = ptrs.ring_ptr->rf_v;
//        ptrs.ring.h = ptrs.ring_ptr->rf_h;
//        ptrs.ring.phi = ptrs.ring_ptr->rf_phi;
//        ptrs.ring.gamma_tr = ptrs.ring_ptr->gamma_tr;

//        ptrs.rf.reset(new RF());
        ptrs.rf->v = ptrs.ring_ptr->rf_v;
        ptrs.rf->h = ptrs.ring_ptr->rf_h;
        ptrs.rf->phi = ptrs.ring_ptr->rf_phi;
        ptrs.rf->gamma_tr = ptrs.ring_ptr->gamma_tr;
//        ptrs.ring->rf = ptrs.rf.get();
    }
    std::cout<<"Ring created!"<<std::endl;
}

void create_cooler(Set_ptrs &ptrs) {
    double length = ptrs.cooler_ptr->length;
    int section_number = ptrs.cooler_ptr->section_number;
    double magnetic_field = ptrs.cooler_ptr->magnetic_field;
    double bet_x = ptrs.cooler_ptr->bet_x;
    double bet_y = ptrs.cooler_ptr->bet_y;
    double disp_x = ptrs.cooler_ptr->disp_x;
    double disp_y = ptrs.cooler_ptr->disp_y;
    double alpha_x = ptrs.cooler_ptr->alpha_x;
    double alpha_y = ptrs.cooler_ptr->alpha_y;
    double disp_dx = ptrs.cooler_ptr->disp_dx;
    double disp_dy = ptrs.cooler_ptr->disp_dy;
    assert(length>0 && section_number>0 && bet_x>0 && bet_y>0 && "WRONG PARAMETER VALUE FOR COOLER!");
    ptrs.cooler.reset(new Cooler(length, section_number, magnetic_field, bet_x, bet_y, disp_x, disp_y, alpha_x, alpha_y,
                                 disp_dx, disp_dy));
    std::cout<<"Cooler created!"<<std::endl;
}

void calculate_ibs(Set_ptrs &ptrs) {
    assert(ptrs.ring.get()!=nullptr && "MUST CREATE THE RING BEFORE CALCULATE IBS RATE!");
    assert(ptrs.ibs_ptr.get()!=nullptr && "PLEASE SET UP THE PARAMETERS FOR IBS RATE CALCULATION!");
    int nu = ptrs.ibs_ptr->nu;
    int nv = ptrs.ibs_ptr->nv;
    int nz = ptrs.ibs_ptr->nz;
    double log_c = ptrs.ibs_ptr->log_c;
    double k = ptrs.ibs_ptr->coupling;
    double rx, ry, rz;

    if (log_c>0) {
        assert(nu>0 && nv>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION!");
        IBSParas ibs_paras(nu, nv, log_c);
        if (k>0) ibs_paras.set_k(k);
        ibs_rate(*ptrs.ring->lattice_, *ptrs.ion_beam, ibs_paras, rx, ry, rz);
    }
    else {
        assert(nu>0 && nv>0 && nz>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION!");
        IBSParas ibs_paras(nu, nv, nz);
        if (k>0) ibs_paras.set_k(k);
        ibs_rate(*ptrs.lattice, *ptrs.ion_beam, ibs_paras, rx, ry, rz);
    }
    ptrs.ibs_rate.at(0) = rx;
    ptrs.ibs_rate.at(1) = ry;
    ptrs.ibs_rate.at(2) = rz;
    std::cout<<std::scientific;
    std::cout << std::setprecision(3);
    std::cout<<"IBS rate (1/s): "<<rx<<"  "<<ry<<"  "<<rz<<std::endl;
}

void calculate_ecool(Set_ptrs &ptrs) {
    assert(ptrs.cooler.get()!=nullptr && "MUST CREATE THE COOLER BEFORE CALCULATE ELECTRON COOLING RATE!");
    assert(ptrs.e_beam.get()!=nullptr && "MUST CREATE THE ELECTRON BEAM BEFORE CALCULATE ELECTRON COOLING RATE!");
    assert(ptrs.ecool_ptr.get()!=nullptr && "PLEASE SET UP THE PARAMETERS FOR ELECTRON COOLING RATE CALCULATION!");
    int n_sample = ptrs.ecool_ptr->n_sample;
    assert(n_sample > 0 && "WRONG PARAMETER VALUE FOR ELECTRON COOLING RATE CALCULATION!");
    EcoolRateParas ecool_paras(n_sample);
    std::string force_formula = ptrs.ecool_ptr->force;
    assert(std::find(FRICTION_FORCE_FORMULA.begin(),FRICTION_FORCE_FORMULA.end(),force_formula)!=FRICTION_FORCE_FORMULA.end()
               && "UNKNOWN FRICTION FORCE FORMULA SECTION_ECOOL!");
    ForceFormula force;
    if (force_formula == "PARKHOMCHUK") {
        force = ForceFormula::PARKHOMCHUK;
    }
    ForceParas force_paras(force);

    assert(ptrs.ion_beam.get()!=nullptr && "MUST CREATE THE ION BEAM BEFORE CALCULATE ELECTRON COOLING RATE!");
    assert(ptrs.ring.get()!=nullptr && "MUST CREATE THE RING BEFORE CALCULATE ELECTRON COOLING RATE!");
    double rx, ry, rz;
    ecooling_rate(ecool_paras, force_paras, *ptrs.ion_beam, *ptrs.cooler, *ptrs.e_beam, *ptrs.ring, rx, ry, rz);
    ptrs.ecool_rate.at(0) = rx;
    ptrs.ecool_rate.at(1) = ry;
    ptrs.ecool_rate.at(2) = rz;
    std::cout<<std::scientific;
    std::cout << std::setprecision(3);
    std::cout<<"Electron cooling rate (1/s): "<<rx<<"  "<<ry<<"  "<<rz<<std::endl;
}

void total_expansion_rate(Set_ptrs &ptrs) {
    if (std::all_of(ptrs.ibs_rate.begin(), ptrs.ibs_rate.end(), [](double i) { return i==0; })) calculate_ibs(ptrs);
    if (std::all_of(ptrs.ecool_rate.begin(), ptrs.ecool_rate.end(), [](double i) { return i==0; })) calculate_ecool(ptrs);
    for(int i=0; i<3; ++i) ptrs.total_rate.at(i) = ptrs.ecool_rate.at(i) + ptrs.ibs_rate.at(i);
    std::cout<<std::scientific;
    std::cout << std::setprecision(3);
    std::cout<<"Total expansion rate (1/s): "<<ptrs.total_rate.at(0)<<"  "<<ptrs.total_rate.at(1)<<"  "
        <<ptrs.total_rate.at(2)<<std::endl;
}

void run_simulation(Set_ptrs &ptrs) {
    assert(ptrs.dynamic_ptr.get()!=nullptr && "PLEASE SET UP THE PARAMETERS FOR SIMULATION!");
    bool ibs = ptrs.dynamic_ptr->ibs;
    bool ecool = ptrs.dynamic_ptr->ecool;
    double t = ptrs.dynamic_ptr->time;
    int n_step = ptrs.dynamic_ptr->n_step;
    int n_sample = ptrs.dynamic_ptr->n_sample;
    int output_intvl = ptrs.dynamic_ptr->output_intvl;
    int save_ptcl_intvl = ptrs.dynamic_ptr->save_ptcl_intvl;

    if(ptrs.dynamic_ptr->model == DynamicModel::TURN_BY_TURN) {
        double dt = ptrs.ring->circ()/(ptrs.ion_beam->beta()*k_c);
        if (n_step>0) {
            t = n_step * dt;
        }
        else if (t>0) {
            n_step = ceil(t/dt);
            t = n_step * dt;
        }
    }

    assert(t>0 && n_step>0 && output_intvl>0 && "WRONG PARAMETERS FOR SIMULAITON!");
    dynamic_paras = new DynamicParas(t, n_step, ibs, ecool);
    dynamic_paras->set_model(ptrs.dynamic_ptr->model);
    dynamic_paras->set_n_sample(n_sample);
    if (output_intvl>1) dynamic_paras->set_output_intvl(output_intvl);
    if (save_ptcl_intvl>0) dynamic_paras->set_ion_save(save_ptcl_intvl);
    dynamic_paras->set_model(ptrs.dynamic_ptr->model);
    dynamic_paras->set_output_file(ptrs.dynamic_ptr->filename);

    if(dynamic_paras->model()==DynamicModel::PARTICLE && !ecool)
        assert(n_sample>0 && "SET N_SAMPLE FOR SIMULATION!");

    assert(ptrs.ring.get()!=nullptr && "MUST CREATE THE RING BEFORE SIMULATION!");
    if (ecool) {
        assert(ptrs.e_beam.get()!=nullptr && "NEED TO CREATE THE E_BEAM BEFORE SIMULATION!");
        assert(ptrs.cooler.get()!=nullptr && "NEED TO CREATE THE COOLER BEFORE SIMULATION!");
        assert(ptrs.ecool_ptr.get()!=nullptr && "PLEASE SET UP THE PARAMETERS FOR ELECTRON COOLING RATE CALCULATION!");
        if (n_sample<=0) {
            n_sample = ptrs.ecool_ptr->n_sample;
            assert(n_sample > 0 && "WRONG PARAMETER VALUE FOR ELECTRON COOLING RATE CALCULATION!");
            dynamic_paras->set_n_sample(n_sample);
        }
        ecool_paras = new EcoolRateParas(n_sample);
        std::string force_formula = ptrs.ecool_ptr->force;
        assert(std::find(FRICTION_FORCE_FORMULA.begin(),FRICTION_FORCE_FORMULA.end(),force_formula)!=FRICTION_FORCE_FORMULA.end()
                   && "UNKNOWN FRICTION FORCE FORMULA SECTION_ECOOL!");
        ForceFormula force;
        if (force_formula == "PARKHOMCHUK") {
            force = ForceFormula::PARKHOMCHUK;
        }
        force_paras = new ForceParas(force);
    }

    if(ibs) {
        assert(ptrs.ibs_ptr.get()!=nullptr && "PLEASE SET UP THE PARAMETERS FOR IBS RATE CALCULATION!");
        int nu = ptrs.ibs_ptr->nu;
        int nv = ptrs.ibs_ptr->nv;
        int nz = ptrs.ibs_ptr->nz;
        double log_c = ptrs.ibs_ptr->log_c;
        double k = ptrs.ibs_ptr->coupling;
        double rx, ry, rz;

        if (log_c>0) {
            assert(nu>0 && nv>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION!");
            ibs_paras = new IBSParas(nu, nv, log_c);
        }
        else {
            assert(nu>0 && nv>0 && nz>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION!");
            ibs_paras = new IBSParas(nu, nv, nz);
        }
        if (k>0) ibs_paras->set_k(k);
    }

    if (ibs && !ecool && dynamic_paras->model()==DynamicModel::PARTICLE) {
        assert(ptrs.dynamic_ptr->ref_bet_x>0 && ptrs.dynamic_ptr->ref_bet_y>0 && "WRONG VALUE FOR REFERENCE TWISS PARAMETERS");
        dynamic_paras->twiss_ref.bet_x = ptrs.dynamic_ptr->ref_bet_x;
        dynamic_paras->twiss_ref.bet_y = ptrs.dynamic_ptr->ref_bet_y;
        dynamic_paras->twiss_ref.alf_x = ptrs.dynamic_ptr->ref_alf_x;
        dynamic_paras->twiss_ref.alf_y = ptrs.dynamic_ptr->ref_alf_y;
        dynamic_paras->twiss_ref.disp_x = ptrs.dynamic_ptr->ref_disp_x;
        dynamic_paras->twiss_ref.disp_y = ptrs.dynamic_ptr->ref_disp_y;
        dynamic_paras->twiss_ref.disp_dx = ptrs.dynamic_ptr->ref_disp_dx;
        dynamic_paras->twiss_ref.disp_dy = ptrs.dynamic_ptr->ref_disp_dy;
    }

    dynamic(*ptrs.ion_beam, *ptrs.cooler, *ptrs.e_beam, *ptrs.ring);
}

void run(std::string &str, Set_ptrs &ptrs) {
    str = trim_blank(str);
    str = trim_tab(str);
    assert(std::find(RUN_COMMANDS.begin(),RUN_COMMANDS.end(),str)!=RUN_COMMANDS.end() && "WRONG COMMANDS IN SECTION_RUN!");
    if (str == "CREATE_ION_BEAM") {
        create_ion_beam(ptrs);
    }
    else if(str == "CREATE_RING") {
        create_ring(ptrs);
    }
    else if(str == "CALCULATE_IBS") {
        calculate_ibs(ptrs);
    }
    else if(str == "CREATE_COOLER") {
        create_cooler(ptrs);
    }
    else if(str == "CREATE_E_BEAM") {
        create_e_beam(ptrs);
    }
    else if(str == "CALCULATE_ECOOL") {
        calculate_ecool(ptrs);
    }
    else if(str == "TOTAL_EXPANSION_RATE") {
        total_expansion_rate(ptrs);
    }
    else if(str == "RUN_SIMULATION") {
        run_simulation(ptrs);
    }
    else {
        assert(false&&"Wrong arguments in section_run!");
    }
}

void set_section_run(Set_ptrs &ptrs) {
    std::fill(ptrs.ibs_rate.begin(), ptrs.ibs_rate.end(), 0);
    std::fill(ptrs.ecool_rate.begin(), ptrs.ecool_rate.end(), 0);
    std::fill(ptrs.total_rate.begin(), ptrs.total_rate.end(), 0);
}

void define_ring(string &str, Set_ring *ring_args) {
    assert(ring_args!=nullptr && "SECTION_RING MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_RING!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    str_toupper(var);
    assert(std::find(RING_ARGS.begin(),RING_ARGS.end(),var)!=RING_ARGS.end() && "WRONG COMMANDS IN SECTION_RING!");
    if (var=="LATTICE") {
        ring_args->lattice_file = val;
    }
    else {
       if(math_parser==NULL) {
            if (var=="QX") {
                ring_args->qx = std::stod(val);
            }
            else if (var=="QY") {
                ring_args->qy = std::stod(val);
            }
            else if (var=="QS") {
                ring_args->qs = std::stod(val);
            }
            else if (var=="GAMMA_TR") {
                ring_args->gamma_tr = std::stod(val);
            }
            else if (var=="RF_V") {
                ring_args->rf_v = std::stod(val);
            }
            else if (var=="RF_H") {
                ring_args->rf_h = std::stoi(val);
            }
            else if (var=="RF_PHI") {
                ring_args->rf_phi = std::stod(val);
            }
            else {
                assert(false&&"Wrong arguments in section_ring!");
            }
       }
       else {
            mupSetExpr(math_parser, val.c_str());
            if (var=="QX") {
                ring_args->qx = mupEval(math_parser);
            }
            else if (var=="QY") {
                ring_args->qy = mupEval(math_parser);
            }
            else if (var=="QS") {
                ring_args->qs = mupEval(math_parser);
            }
            else if (var=="GAMMA_TR") {
                ring_args->gamma_tr = mupEval(math_parser);
            }
            else if (var=="RF_V") {
                ring_args->rf_v = mupEval(math_parser);
            }
            else if (var=="RF_H") {
                ring_args->rf_h = mupEval(math_parser);
            }
            else if (var=="RF_PHI") {
                ring_args->rf_phi = mupEval(math_parser);
            }
            else {
                assert(false&&"Wrong arguments in section_ring!");
            }
       }

    }
}

void define_cooler(std::string &str, Set_cooler *cooler_args) {
    assert(cooler_args!=nullptr && "SECTION_COOLER MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_RING!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(COOLER_ARGS.begin(),COOLER_ARGS.end(),var)!=COOLER_ARGS.end() && "WRONG COMMANDS IN SECTION_RING!");
    if(math_parser==NULL) {
        if (var=="LENGTH") {
            cooler_args->length = std::stod(val);
        }
        else if (var == "SECTION_NUMBER") {
            cooler_args->section_number = std::stoi(val);
        }
        else if (var == "MAGNETIC_FIELD") {
            cooler_args->magnetic_field = std::stod(val);
        }
        else if (var == "BET_X") {
            cooler_args->bet_x = std::stod(val);
        }
        else if (var == "BET_Y") {
            cooler_args->bet_y = std::stod(val);
        }
        else if (var == "DISP_X") {
            cooler_args->disp_x = std::stod(val);
        }
        else if (var == "DISP_Y") {
            cooler_args->disp_y = std::stod(val);
        }
        else if (var == "ALPHA_X") {
            cooler_args->alpha_x = std::stod(val);
        }
        else if (var == "ALPHA_Y") {
            cooler_args->alpha_y = std::stod(val);
        }
        else if (var == "DISP_DX") {
            cooler_args->disp_dx = std::stod(val);
        }
        else if (var == "DISP_DY") {
            cooler_args->disp_dy = std::stod(val);
        }
        else {
            assert(false&&"Wrong arguments in section_cooler!");
        }
    }
    else {
        mupSetExpr(math_parser, val.c_str());
        if (var=="LENGTH") {
            cooler_args->length = mupEval(math_parser);
        }
        else if (var == "SECTION_NUMBER") {
            cooler_args->section_number = static_cast<int>(mupEval(math_parser));
        }
        else if (var == "MAGNETIC_FIELD") {
            cooler_args->magnetic_field = mupEval(math_parser);
        }
        else if (var == "BET_X") {
            cooler_args->bet_x = mupEval(math_parser);
        }
        else if (var == "BET_Y") {
            cooler_args->bet_y = mupEval(math_parser);
        }
        else if (var == "DISP_X") {
            cooler_args->disp_x = mupEval(math_parser);
        }
        else if (var == "DISP_Y") {
            cooler_args->disp_y = mupEval(math_parser);
        }
        else if (var == "ALPHA_X") {
            cooler_args->alpha_x = mupEval(math_parser);
        }
        else if (var == "ALPHA_Y") {
            cooler_args->alpha_y = mupEval(math_parser);
        }
        else if (var == "DISP_DX") {
            cooler_args->disp_dx = mupEval(math_parser);
        }
        else if (var == "DISP_DY") {
            cooler_args->disp_dy = mupEval(math_parser);
        }
        else {
            assert(false&&"Wrong arguments in section_cooler!");
        }
    }
}

void set_ibs(string &str, Set_ibs *ibs_args) {
    assert(ibs_args!=nullptr && "SECTION_IBS MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_IBS!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(IBS_ARGS.begin(),IBS_ARGS.end(),var)!=IBS_ARGS.end() && "WRONG COMMANDS IN SECTION_IBS!");
    if (math_parser == NULL) {
        if (var == "NU") {
            ibs_args->nu = std::stoi(val);
        }
        else if(var == "NV") {
            ibs_args->nv = std::stoi(val);
        }
        else if(var == "NZ") {
            ibs_args->nz = std::stoi(val);
        }
        else if(var == "LOG_C") {
            ibs_args->log_c = std::stod(val);
        }
        else if(var == "COUPLING") {
            ibs_args->coupling = std::stod(val);
        }
        else {
            assert(false&&"Wrong arguments in section_ibs!");
        }
    }
    else {
        mupSetExpr(math_parser, val.c_str());
        if (var == "NU") {
            ibs_args->nu = static_cast<int>(mupEval(math_parser));
        }
        else if(var == "NV") {
            ibs_args->nv = static_cast<int>(mupEval(math_parser));
        }
        else if(var == "NZ") {
            ibs_args->nz = static_cast<int>(mupEval(math_parser));
        }
        else if(var == "LOG_C") {
            ibs_args->log_c = mupEval(math_parser);
        }
        else if(var == "COUPLING") {
            ibs_args->coupling = mupEval(math_parser);
        }
        else {
            assert(false&&"Wrong arguments in section_ibs!");
        }
    }

}

void set_simulation(string &str, Set_dynamic *dynamic_args) {
    assert(dynamic_args!=nullptr && "SECTION_SIMULATION MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_SIMULATION!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(SIMULATION_ARGS.begin(),SIMULATION_ARGS.end(),var)!=SIMULATION_ARGS.end() && "WRONG COMMANDS IN SECTION_SIMULATION!");

    if (var == "MODEL") {
        if (val == "RMS") dynamic_args->model = DynamicModel::RMS;
        else if(val == "MODEL_BEAM") dynamic_args->model = DynamicModel::MODEL_BEAM;
        else if(val == "PARTICLE") dynamic_args->model = DynamicModel::PARTICLE;
        else if(val == "TURN_BY_TURN") dynamic_args->model = DynamicModel::TURN_BY_TURN;
        else assert("DYNAMIC MODEL NOT SUPPORTED IN SIMULATION!");
    }
    else if (var == "OUTPUT_FILE") {
        dynamic_args->filename = val;
    }
    else if (var == "IBS" ) {
        if (val == "ON") dynamic_args->ibs = true;
        else if (val == "OFF") dynamic_args->ibs = false;
        else assert("WRONG VALUE IN SECTION_SIMULATION!");
    }
    else if (var == "E_COOL") {
        if (val == "ON") dynamic_args->ecool = true;
        else if (val == "OFF") dynamic_args->ecool = false;
        else assert("WRONG VALUE IN SECTION_SIMULATION!");
    }
    else {
        if (math_parser == NULL) {
            if (var == "TIME") {
                dynamic_args->time = std::stod(val);
            }
            else if (var == "STEP_NUMBER") {
                dynamic_args->n_step = std::stoi(val);
            }
            else if (var == "SAMPLE_NUMBER") {
                dynamic_args->n_sample = std::stoi(val);
            }
            else if (var == "OUTPUT_INTERVAL") {
                dynamic_args->output_intvl = std::stoi(val);
            }
            else if (var == "SAVE_PARTICLE_INTERVAL") {
                dynamic_args->save_ptcl_intvl = std::stoi(val);
            }
            else if (var == "REF_BET_X") {
                dynamic_args->ref_bet_x = std::stod(val);
            }
            else if (var == "REF_BET_Y") {
                dynamic_args->ref_bet_y = std::stod(val);
            }
            else if (var == "REF_ALF_X") {
                dynamic_args->ref_alf_x = std::stod(val);
            }
            else if (var == "REF_ALF_Y") {
                dynamic_args->ref_alf_y = std::stod(val);
            }
            else if (var == "REF_DISP_X") {
                dynamic_args->ref_disp_x = std::stod(val);
            }
            else if (var == "REF_DISP_Y") {
                dynamic_args->ref_disp_y = std::stod(val);
            }
            else if (var == "REF_DISP_DX") {
                dynamic_args->ref_disp_dx = std::stod(val);
            }
            else if (var == "REF_DISP_DY") {
                dynamic_args->ref_disp_dy = std::stod(val);
            }
            else {
                assert(false&&"Wrong arguments in section_simulation!");
            }
        }
        else {
            mupSetExpr(math_parser, val.c_str());
            if (var == "TIME") {
                dynamic_args->time = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "STEP_NUMBER") {
                dynamic_args->n_step = static_cast<int>(mupEval(math_parser));
            }
            else if (var == "SAMPLE_NUMBER") {
                dynamic_args->n_sample = static_cast<int>(mupEval(math_parser));
            }
            else if (var == "OUTPUT_INTERVAL") {
                dynamic_args->output_intvl = static_cast<int>(mupEval(math_parser));
            }
            else if (var == "SAVE_PARTICLE_INTERVAL") {
                dynamic_args->save_ptcl_intvl = static_cast<int>(mupEval(math_parser));
            }
            else if (var == "REF_BET_X") {
                dynamic_args->ref_bet_x = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_BET_Y") {
                dynamic_args->ref_bet_y = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_ALF_X") {
                dynamic_args->ref_alf_x = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_ALF_Y") {
                dynamic_args->ref_alf_y = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_DISP_X") {
                dynamic_args->ref_disp_x = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_DISP_Y") {
                dynamic_args->ref_disp_y = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_DISP_DX") {
                dynamic_args->ref_disp_dx = static_cast<double>(mupEval(math_parser));
            }
            else if (var == "REF_DISP_DY") {
                dynamic_args->ref_disp_dy = static_cast<double>(mupEval(math_parser));
            }
            else {
                assert(false&&"Wrong arguments in section_simulation!");
            }
        }
    }
}

void set_ecool(string &str, Set_ecool *ecool_args){
    assert(ecool_args!=nullptr && "SECTION_ECOOL MUST BE CLAIMED!");
    string::size_type idx = str.find("=");
    assert(idx!=string::npos && "WRONG COMMAND IN SECTION_ECOOL!");
    string var = str.substr(0, idx);
    string val = str.substr(idx+1);
    var = trim_blank(var);
    var = trim_tab(var);
    val = trim_blank(val);
    val = trim_tab(val);
    assert(std::find(ECOOL_ARGS.begin(),ECOOL_ARGS.end(),var)!=ECOOL_ARGS.end() && "WRONG COMMANDS IN SECTION_ECOOL!");

    if (var == "FORCE_FORMULA") {
        ecool_args->force = val;
    }
    else {
        if (math_parser == NULL) {
            if (var == "SAMPLE_NUMBE") {
                ecool_args->n_sample = std::stod(val);
            }
            else {
                assert(false&&"Wrong arguments in section_ecool!");
            }
        }
        else {
            mupSetExpr(math_parser, val.c_str());
            if (var == "SAMPLE_NUMBER") {
                ecool_args->n_sample = static_cast<double>(mupEval(math_parser));
            }
            else {
                assert(false&&"Wrong arguments in section_ecool!");
            }
        }
    }
}

void parse(std::string &str, muParserHandle_t &math_parser){

    if (str == "LIST_VAR") {
        ListVar(math_parser);
    }
    else if(str == "LIST_CONST") {
        ListConst(math_parser);
    }
    else if (str.substr(0,5) == "PRINT") {
        string var = str.substr(6);
        var = trim_blank(var);
        var = trim_tab(var);
        mupSetExpr(math_parser, var.c_str());
        std::cout<<var<<" = "<<mupEval(math_parser)<<std::endl;
    }
    else {
        mupSetExpr(math_parser, str.c_str());
        mupEval(math_parser);
    }
}
