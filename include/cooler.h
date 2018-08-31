#ifndef COOLER_H
#define COOLER_H

class Cooler{
    double length_;      // in meter
    double section_number_;
    double magnetic_field_;           // in Tesla
    double beta_h_;      // in meter
    double beta_v_;      // in meter
    double disp_h_;
    double disp_v_;
    double alpha_h_;
    double alpha_v_;
    double der_disp_h_;
    double der_disp_v_;
 public:
    double length(){return length_;}
    double section_number(){return section_number_;}
    double magnetic_field(){return magnetic_field_;}
    double beta_h(){return beta_h_;}
    double beta_v(){return beta_v_;}
    double alpha_h(){return alpha_h_;}
    double alpha_v(){return alpha_v_;}
    double disp_h(){return disp_h_;}
    double disp_v(){return disp_v_;}
    double der_disp_h(){return der_disp_h_;}
    double der_disp_v(){return der_disp_v_;}
    Cooler(double length, double section_number, double magnetic_field, double beta_h, double beta_v, double disp_h=0,
           double disp_v=0, double alpha_h=0, double alpha_v=0, double der_disp_h=0, double der_disp_v=0):length_(length),
           section_number_(section_number), magnetic_field_(magnetic_field), beta_h_(beta_h), beta_v_(beta_v),
           disp_h_(disp_h), disp_v_(disp_v), alpha_h_(alpha_h), alpha_v_(alpha_v), der_disp_h_(der_disp_h),
           der_disp_v_(der_disp_v){};
};

#endif // COOLER_H
