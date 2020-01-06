#ifndef LIUMINOSITY_H
#define LIUMINOSITY_H

#include <vector>

struct CollidingBeam {
    double np = 0;
    double sigma_x = 0;
    double sigma_y = 0;
    double geo_emit_x = 0;
    double geo_emit_y = 0;
    double bet_x_star = 0;
    double bet_y_star = 0;
    double bet_x_max = 0;
    double bet_y_max = 0;
    double focus_length_x = 0;
    double focus_length_y = 0;
    double aper_x = 0;
    double aper_y = 0;
    bool adjust_bet = false;
};

class Luminosity {
    std::vector<CollidingBeam> beam = std::vector<CollidingBeam>(2);
    double dx_ = 0;
    double dy_ = 0;
    double freq_ = 1;
    double aper_ratio_ = 10;
    bool match_= false;
    bool use_ion_emit_ = false;
public:
    void set_distance(double dx, double dy){dx_=dx; dy_=dy;}
    void set_freq(double f){freq_=f;}
    void set_match(bool b){match_ = b;}
    void set_geo_emit(double emit_x, double emit_y, int i);
    void set_beam_size(double sigma_x, double sigma_y, int i);
    void set_particle_number(double n, int i);
    void set_bet(double bet_x, double bet_y, int i);
    bool match(){return match_;}
    void set_focus_length(double x, double y, int i);
    void set_aperture(double x, double y, int i);
    void calc_bet_max(int i);
    void adjust_bet(int i);
    void set_adjust_bet(bool b, int i);
    void set_aper_ratio(double i){aper_ratio_ = i;}
    void set_use_ion_emit(bool b){use_ion_emit_ = b;}
    bool use_ion_emit(){return use_ion_emit_;}
    void match(int i, int j); //match the size of beam j to that of beam i.
    double bet_x_star(int i);
    double bet_y_star(int i);
    double bet_x_max(int i);
    double bet_y_max(int i);
    double luminosity();
};
//
//class Luminosity {
//    double dx_ = 0;
//    double dy_ = 0;
//    double np_1_ = 0;
//    double np_2_ = 0;
//    double freq_ = 1;
//    double sigma_x1_ = 0;
//    double sigma_y1_ = 0;
//    double sigma_x2_ = 0;
//    double sigma_y2_ = 0;
//    double bet_x1_ = 0;
//    double bet_y1_ = 0;
//    double bet_x2_ = 0;
//    double bet_y2_ = 0;
//    double geo_emit_x1_ = 0;
//    double geo_emit_y1_ = 0;
//    double geo_emit_x2_ = 0;
//    double geo_emit_y2_ = 0;
//    bool use_ion_emittance_ = true;
//    double bet_x1_max_ = 0;
//    double bet_y1_max_ = 0;
//    double focus_length_x1_ = 0;
//    double focus_length_y1_ = 0;
//    double aper_x1_ = 0;
//    double aper_y1_ = 0;
//    double bet_x2_max_ = 0;
//    double bet_y2_max_ = 0;
//    double focus_length_x2_ = 0;
//    double focus_length_y2_ = 0;
//    double aper_x2_ = 0;
//    double aper_y2_ = 0;
//    bool adjust_bet_ = false;
//public:
//    void set_distance(double dx, double dy){dx_=dx; dy_=dy;}
//    void set_freq(double f){freq_=f;}
//    void set_use_ion_emit(bool b){use_ion_emittance_ = b;}
//    void set_geo_emit(double emit_x, double emit_y, int i);
//    void set_beam_size(double sigma_x, double sigma_y, int i);
//    void set_particle_number(double n, int i);
//    void set_bet(double bet_x, double bet_y, int i);
//    bool use_ion_emittance(){return use_ion_emittance_;}
//    void set_focus_length(double x, double y, int i);
//    void set_aperture(double x, double y, int i);
//    void calc_bet_max(int i);
//    void adjust_bet(int i);
//    void set_adjust_bet(bool b) {adjust_bet_=b;}
//    bool adjust_bet(){return adjust_bet_;}
//
//    double luminosity();
//};

#endif // LIUMINOSITY_H
