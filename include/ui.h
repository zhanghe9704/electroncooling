#ifndef UI_H_INCLUDED
#define UI_H_INCLUDED

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "beam.h"
#include "cooler.h"
#include "dynamic.h"
#include "ecooling.h"
#include "math_parser.h"
#include "ring.h"

//enum class sections{SECTION_ION, SECTION_RING, SECTION_COOLER};
//enum class section_ion{CHARGE_NUMBER, MASS, KINETIC_ENERGY, NORM_EMIT_X, NORM_EMIT_Y, MOMENTUM_SPREAD, PARTICLE_NUMBER,
//    RMS_BUNCH_LENGHT, DEFINE_ION_BEAM};



struct Set_ion{
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

struct Set_ring{
    std::string lattice_file = "";
    double qx = 0;
    double qy = 0;
    double qs = 0;
    double rf_v = 0;
    double rf_phi = 0;
    int rf_h = 1;
    double gamma_tr = 0;
};

struct Set_ibs{
     int nu = 0;
     int nv = 0;
     int nz = 0;
     double log_c = 0;
     double coupling = -1;
     IBSModel model = IBSModel::MARTINI;
     bool ibs_by_element = false;
};

struct Set_ecool{
     int n_sample = 0;
//     std::string force = "PARKHOMCHUK";
     ForceFormula force = ForceFormula::PARKHOMCHUK;
};

struct Set_cooler{
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


struct Set_e_beam{
     double gamma = 0;
     double tmp_tr = 0;
     double tmp_l = 0;
     double n = 0;
     double sigma_x = 0;
     double sigma_y = 0;
     double sigma_z = 0;
     double current = 0;
     double radius = 0;
     double length = 0;
     double rh = 0;
     double rv = 0;
     double r_inner = 0;
     double r_outter = 0;
     int line_skip = 0;
     int n_particle = 0;
     int particle_perbox = 200;
     std::string shape = "";
     std::string particle_file = "";
     bool corr = false;
     bool binary = false;
     int buffer = 1000;
};

struct Set_dynamic{
    double time = 0;
    int n_step = 0;
    int n_sample = 0;
    bool ibs = true;
    bool ecool = true;
    bool fixed_bunch_length = false;
    bool reset_time = true;
    bool overwrite = true;
    bool calc_luminosity = false;
    int output_intvl = 1;
    int save_ptcl_intvl = -1;
    std::string filename = "output_dynamic.txt";
    DynamicModel model = DynamicModel::RMS;
    double ref_bet_x = 0;
    double ref_bet_y = 0;
    double ref_alf_x = 0;
    double ref_alf_y = 0;
    double ref_disp_x = 0;
    double ref_disp_y = 0;
    double ref_disp_dx = 0;
    double ref_disp_dy = 0;
};

struct Set_luminosity{
    double dx = 0;
    double dy = 0;
    double np_1 = 0;
    double np_2 = 0;
    double freq = 1;
    double sigma_x1 = 0;
    double sigma_y1 = 0;
    double sigma_x2 = 0;
    double sigma_y2 = 0;
    double bet_x1 = 0;
    double bet_y1 = 0;
    double bet_x2 = 0;
    double bet_y2 = 0;
    double geo_emit_x1 = 0;
    double geo_emit_x2 = 0;
    double geo_emit_y1 = 0;
    double geo_emit_y2 = 0;
    bool use_ion_emittance = true;
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
     std::unique_ptr<Set_dynamic> dynamic_ptr = nullptr;
     std::unique_ptr<Set_luminosity> luminosity_ptr = nullptr;
//     std::unique_ptr<Tunes> tunes = nullptr;
//     std::unique_ptr<RF> rf = nullptr;
};

enum class Section{NONE, SECTION_ION, SECTION_RING, SECTION_COOLER, SECTION_RUN, SECTION_IBS, SECTION_SCRATCH,
    SECTION_E_BEAM_SHAPE, SECTION_E_BEAM, SECTION_ECOOL, SECTION_SIMULATION, SECTION_LUMINOSITY};

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
void set_luminosity(string &str, Set_luminosity *lum_args);
void set_ecool(std::string &str, Set_ecool *ecool_args);
void set_section_run(Set_ptrs &ptrs);
void set_simulation(std::string &str, Set_dynamic *dynamic_args);
#endif // UI_H_INCLUDED
