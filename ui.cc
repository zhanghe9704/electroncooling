#include "ui.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "beam.h"
#include "constants.h"
#include "ibs.h"

using std::string;
muParserHandle_t math_parser = NULL;
//std::vector<std::string> sections = {"SECTION_ION", "SECTION_RING", "SECTION_COOLER"};
std::vector<string> ION_ARGS = {"CHARGE_NUMBER", "MASS", "KINETIC_ENERGY", "NORM_EMIT_X", "NORM_EMIT_Y",
    "MOMENTUM_SPREAD", "PARTICLE_NUMBER", "RMS_BUNCH_LENGTH"};
std::vector<string> RUN_COMMANDS = {"CREATE_ION_BEAM", "CREATE_RING", "CREATE_ELEC_BEAM", "CREATE_COOLER",
    "CALCULATE_IBS"};
std::vector<string> RING_ARGS = {"LATTICE"};
std::vector<string> IBS_ARGS = {"NU","NV","NZ","LOG_C","COUPLING"};
std::vector<string> COOLER_ARGS = {"LENGTH", "SECTION_NUMBER", "MAGNETIC_FIELD", "BET_X", "BET_Y", "DISP_X", "DISP_Y",
    "ALPHA_X", "ALPHA_Y", "DISP_DX", "DISP_DY"};
//std::vector<string> SCRATCH_COMMANDS = {"PRINT", "LIST_VAR", "LIST_CONST"};

std::map<std::string, Section> sections{
    {"SECTION_ION",Section::SECTION_ION},
    {"SECTION_RING",Section::SECTION_RING},
    {"SECTION_COOLER",Section::SECTION_COOLER},
    {"SECTION_RUN",Section::SECTION_RUN},
    {"SECTION_IBS",Section::SECTION_IBS},
    {"SECTION_SCRATCH", Section::SECTION_SCRATCH}
    };

//enum class Section {
//    SECTION_ION = sections["SECTION_ION"],
//    SECTION_RING = sections["SECTION_RING"],
//    SECTION_COOLER = sections["SECTION_COOLER"]
//    };

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


//void define_ion_beam(std::string &str, Set_ion *ion_args) {
//    assert(ion_args!=nullptr && "SECTION_ION MUST BE CLAIMED!");
//    string::size_type idx = str.find("=");
//    if (idx==string::npos) {
//        //Command
//        str = trim_blank(str);
//        str = trim_tab(str);
//        assert((str=="DEFINE_ION_BEAM") && "Wrong command!");
//    }
//    else {
//        //Set variables
//        string var = str.substr(0, idx);
//        string val = str.substr(idx+1);
//        var = trim_blank(var);
//        var = trim_tab(var);
//        val = trim_blank(val);
//        val = trim_tab(val);
//        std::cout<<var<<" "<<val<<std::endl;
////        assert(std::find(ION_ARGS.begin(), ION_ARGS.end(), var)!=ION_ARGS.end() && "WRONG ARGUMENTS FOR ION BEAM!");
//
//    }
//}

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
    assert(ptrs.ion_beam.get()!=nullptr && "MUST DEFINE THE ION BEFORE DEFINE RING!");
    ptrs.ring.reset(new Ring(*ptrs.lattice, *ptrs.ion_beam));
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
    int nu = ptrs.ibs_ptr->nu;
    int nv = ptrs.ibs_ptr->nv;
    int nz = ptrs.ibs_ptr->nz;
    double log_c = ptrs.ibs_ptr->log_c;
    double rx, ry, rz;

    if (log_c>0) {
        assert(nu>0 && nv>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION");
        IBSParas ibs_paras(nu, nv, log_c);
        ibs_rate(*ptrs.ring->lattice_, *ptrs.ion_beam, ibs_paras, rx, ry, rz);
    }
    else {
        assert(nu>0 && nv>0 && nz>0 && "WRONG PARAMETER VALUE FOR IBS RATE CALCULATION");
        IBSParas ibs_paras(nu, nv, nz);
        ibs_rate(*ptrs.lattice, *ptrs.ion_beam, ibs_paras, rx, ry, rz);
    }
    std::cout<<std::scientific;
    std::cout << std::setprecision(3);
    std::cout<<"IBS rate (1/s): "<<rx<<"  "<<ry<<"  "<<rz<<std::endl;
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
    assert(std::find(RING_ARGS.begin(),RING_ARGS.end(),var)!=RING_ARGS.end() && "WRONG COMMANDS IN SECTION_RING!");
    if (var=="LATTICE") {
        ring_args->lattice_file = val;
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
    }
}

void set_ibs(string &str, Set_ibs *ibs_args) {
    assert(ibs_args!=nullptr && "SECTION_ibs MUST BE CLAIMED!");
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

