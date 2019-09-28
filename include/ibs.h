#ifndef IBS_HPP
#define IBS_HPP

#include <assert.h>
#include <memory>

class Lattice;
class Beam;

enum class IBSModel {MARTINI, BM};

class IBSSolver {
protected:
    int nu_ = 0;                //Grid number in u direction.
    int nv_ = 0;                //Grid number in v direction.
    int nz_ = 0;                //Grid number in z direction.
    double log_c_ = 0.0;     //Coulomb logarithm.
    double k_ = 0.0;          //Coupling rate in transverse directions.
    bool reset_ = true;
    
    void ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y);
public:
    int nu() const { return nu_; }
    int nv() const { return nv_; }
    int nz() const { return nz_; }
    double log_c() const { return log_c_; }
    double k() const { return k_; }
    bool reset() const { return reset_; }
    void set_k(double x) { k_ = x; }
    void set_log_c(double x) { log_c_ = x; }
    void reset_off() { reset_ = false; }
    void reset_on() { reset_ = true; }
    void set_nu(int nu) { assert(nu>0&&"Wrong value of nu in IBS parameters!"); nu_ = nu; }
    void set_nv(int nv) { assert(nv>0&&"Wrong value of nv in IBS parameters!"); nv_ = nv; }
    void set_nz(int nz) { assert(nz>0&&"Wrong value of nz in IBS parameters!"); nz_ = nz; }

    IBSSolver(int nu, int nv, int nz, double log_c, double k);
    
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs) = 0;
};

class IBSSolver_Martini : public IBSSolver {
private:
    // Scratch variables for IBS calculation (Martini model)
    std::unique_ptr<double []> sigma_xbet, sigma_xbetp, sigma_y, sigma_yp;
    std::unique_ptr<double []> a_f, b2_f, c2_f, d2_f, dtld_f, k1_f, k2_f, k3_f;

    std::unique_ptr<double []> sin_u, sin_u2, cos_u2, sin_v, cos_v, sin_u2_cos_v2;
    std::unique_ptr<double []> g1, g2_1, g2_2, g3;
    std::unique_ptr<double []> f1, f2, f3;
    
    void set_bunch_size(int n);
    void bunch_size(const Lattice &lattice, const Beam &beam);
    void set_abcdk(int n);
    void abcdk(const Lattice &lattice, const Beam &beam);
    void coef_f(int nu, int nv);
    void f(int n_element);
    double coef_a(const Lattice &lattice, const Beam &beam) const;
public:
    IBSSolver_Martini(int nu, int nv, int nz, double log_c, double k);
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
};

// IBS calculation by Bjorken-Mtingwa model using Sergei Nagaitsev's method
class IBSSolver_BM : public IBSSolver {
private:
    std::unique_ptr<double []> phi;
    std::unique_ptr<double []> dx2; //D_x * D_x
    std::unique_ptr<double []> dx_betax_phi_2; // D_x * D_x / (beta_x * beta_x) + phi * phi
    std::unique_ptr<double []> sqrt_betay;  // sqrt(beta_y)
    std::unique_ptr<double []> gamma_phi_2; // gamma * gamma * phi * phi
    std::unique_ptr<double []> psi;
    std::unique_ptr<double []> sx, sp, sxp;
    std::unique_ptr<double []> inv_sigma;

    void alloc_var(std::unique_ptr<double []>& ptr, int n);
    void init_fixed_var(const Lattice &lattice, const Beam &beam);
    void calc_kernels(const Lattice &lattice, const Beam &beam);
    double coef_bm(const Lattice &lattice, const Beam &beam) const;
public:
    IBSSolver_BM(double log_c, double k);
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
};

#endif
