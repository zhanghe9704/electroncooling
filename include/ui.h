#ifndef UI_H_INCLUDED
#define UI_H_INCLUDED

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "beam.h"
#include "cooler.h"
#include "ecooling.h"
#include "math_parser.h"
#include "ring.h"

//enum class sections{SECTION_ION, SECTION_RING, SECTION_COOLER};
//enum class section_ion{CHARGE_NUMBER, MASS, KINETIC_ENERGY, NORM_EMIT_X, NORM_EMIT_Y, MOMENTUM_SPREAD, PARTICLE_NUMBER,
//    RMS_BUNCH_LENGHT, DEFINE_ION_BEAM};



class Set_ion{
 public:
    int n_charge = 0;
    double mass = 0;
    double k_energy = 0;
    double emit_nx = 0;
    double emit_ny = 0;
    double dp_p = -1;
    double n_ptcl = 0;
    double ds = 0;
//    bool define_ion_beam = false;
};

class Set_ring{
 public:
    std::string lattice_file = "";
};

class Set_ibs{
 public:
     int nu = 0;
     int nv = 0;
     int nz = 0;
     double log_c = 0;
     double coupling = -1;
};

class Set_ecool{
 public:
     int n_sample = 0;
     std::string force = "PARKHOMCHUK";
};

class Set_cooler{
 public:
     double length = 0;
     int section_number = 0;
     double magnetic_field = 0;
     double bet_x = 0;
     double bet_y = 0;
     double disp_x = 0;
     double disp_y = 0;
     double disp_dx = 0;
     double disp_dy = 0;
     double alpha_x = 0;
     double alpha_y = 0;
};


class Set_e_beam{
 public:
     double gamma = 0;
     double tmp_tr = -1;
     double tmp_l = -1;
     int n = 0;
     double sigma_x = 0;
     double sigma_y = 0;
     double sigma_z = 0;
     double current = 0;
     double radius = 0;
     double length = 0;
     std::string shape = "";
};

class Set_ptrs{
 public:
     std::vector<double> ibs_rate = {0,0,0};
     std::vector<double> ecool_rate = {0,0,0};
     std::vector<double> total_rate = {0,0,0};
     std::unique_ptr<Set_ion> ion_ptr = nullptr;
     std::unique_ptr<Beam> ion_beam = nullptr;
     std::unique_ptr<Set_ring> ring_ptr = nullptr;
     std::unique_ptr<Lattice> lattice = nullptr;
     std::unique_ptr<Set_ibs> ibs_ptr = nullptr;
     std::unique_ptr<Ring> ring = nullptr;
     std::unique_ptr<Set_cooler> cooler_ptr = nullptr;
     std::unique_ptr<Cooler> cooler = nullptr;
     std::unique_ptr<EBeamShape> e_beam_shape = nullptr;
     std::unique_ptr<Set_e_beam> e_beam_ptr = nullptr;
     std::unique_ptr<EBeam> e_beam = nullptr;
     std::unique_ptr<Set_ecool> ecool_ptr = nullptr;
     std::unique_ptr<EcoolRateParas> ecool_paras = nullptr;
     std::unique_ptr<ForceParas> force_paras = nullptr;

};

enum class Section{NONE, SECTION_ION, SECTION_RING, SECTION_COOLER, SECTION_RUN, SECTION_IBS, SECTION_SCRATCH,
    SECTION_E_BEAM_SHAPE, SECTION_E_BEAM, SECTION_ECOOL};

std::string remove_comments(std::string input_line);
std::string trim_blank(std::string input_line);
std::string trim_tab(std::string input_line);
void str_toupper(std::string &str);
void define_ion_beam(std::string &str, Set_ion *ion_args);
void run(std::string &str, Set_ptrs &ptrs);
void define_ring(std::string &str, Set_ring *ring_args);
void set_ibs(std::string &str, Set_ibs *ibs_args);
void parse(std::string &str, muParserHandle_t &math_parser);
void define_cooler(std::string &str, Set_cooler *cooler_args);
void create_cooler(Set_ptrs &ptrs);
void define_e_beam(std::string &str, Set_e_beam *e_beam_args);
void set_ecool(std::string &str, Set_ecool *ecool_args);
#endif // UI_H_INCLUDED
