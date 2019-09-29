#include <fstream>
#include <memory>
#include <cmath>
#include <cstring>
#include "functions.h"
#include "ibs.h"
#include "ring.h"
#include "beam.h"

IBSSolver::IBSSolver(double log_c, double k)
    : log_c_(log_c), k_(k)
{
}

void IBSSolver::ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y)
{
    double rxc = 0.5*(rx*(2-k)+ry*k*emit_y/emit_x);
    double ryc = 0.5*(ry*(2-k)+rx*k*emit_x/emit_y);
    rx = rxc;
    ry = ryc;
}


IBSSolver_Martini::IBSSolver_Martini(int nu, int nv, int nz, double log_c, double k)
    : IBSSolver(log_c, k), nu_(nu), nv_(nv), nz_(nz)
{
#ifndef NDEBUG
	std::cerr << "DEBUG: ISBSolver_Martini constructor" << std::endl;
#endif
}

//Set ptrs for bunch_size
void IBSSolver_Martini::set_bunch_size(int n) {
    if(sigma_xbet.get()==nullptr) {
        sigma_xbet.reset(new double[n]);
        memset(sigma_xbet.get(), 0, n*sizeof(double));
//        std::cout<<"reset bunch_sizze"<<std::endl;
    }
    if(sigma_xbetp.get()==nullptr) {
        sigma_xbetp.reset(new double[n]);
        memset(sigma_xbetp.get(), 0, n*sizeof(double));
    }
    if(sigma_y.get()==nullptr) {
        sigma_y.reset(new double[n]);
        memset(sigma_y.get(), 0, n*sizeof(double));
    }
    if(sigma_yp.get()==nullptr) {
        sigma_yp.reset(new double[n]);
        memset(sigma_yp.get(), 0, n*sizeof(double));
    }
}

//Calculate sigma_xbet, sigma_xbetp, sigma_y, sigma_yp
void IBSSolver_Martini::bunch_size(const Lattice &lattice, const Beam &beam) {
    double emit_x = beam.emit_x();
    double emit_y = beam.emit_y();
    set_bunch_size(lattice.n_element());
    for(int i=0; i<lattice.n_element(); ++i) {
        sigma_xbet[i] = sqrt(lattice.betx(i)*emit_x);
        sigma_y[i] = sqrt(lattice.bety(i)*emit_y);
        double alf2 = lattice.alfx(i);
        alf2 *= alf2;
        sigma_xbetp[i] = sqrt((1+alf2)*emit_x/lattice.betx(i));
        alf2 = lattice.alfy(i);
        alf2 *= alf2;
        sigma_yp[i] = sqrt((1+alf2)*emit_y/lattice.bety(i));
    }
}

//Set the ptrs for abcdk
void IBSSolver_Martini::set_abcdk(int n) {
    if(a_f.get()==nullptr) {
        a_f.reset(new double[n]);
        memset(a_f.get(), 0, n*sizeof(double));
//        std::cout<<"reset abcdk"<<std::endl;
    }
    if(b2_f.get()==nullptr) {
        b2_f.reset(new double[n]);
        memset(b2_f.get(), 0, n*sizeof(double));
    }
    if(c2_f.get()==nullptr) {
        c2_f.reset(new double[n]);
        memset(c2_f.get(), 0, n*sizeof(double));
    }
    if(d2_f.get()==nullptr) {
        d2_f.reset(new double[n]);
        memset(d2_f.get(), 0, n*sizeof(double));
    }
    if(dtld_f.get()==nullptr) {
        dtld_f.reset(new double[n]);
        memset(dtld_f.get(), 0, n*sizeof(double));
    }
    if(k1_f.get()==nullptr) {
        k1_f.reset(new double[n]);
        memset(k1_f.get(), 0, n*sizeof(double));
    }
    if(k2_f.get()==nullptr) {
        k2_f.reset(new double[n]);
        memset(k2_f.get(), 0, n*sizeof(double));
    }
    if(k3_f.get()==nullptr) {
        k3_f.reset(new double[n]);
        memset(k3_f.get(), 0, n*sizeof(double));
    }
}

//Calculate a, b, c, d and dtld
//Call bunch_size() before calling this one
void IBSSolver_Martini::abcdk(const Lattice &lattice, const Beam &beam)
{
    double d_tld, q, sigma_x, sigma_tmp;
    int n = lattice.n_element();

    set_abcdk(n);

    double dp_p = beam.dp_p();
    double beta = beam.beta();
    double gamma = beam.gamma();
    double r = beam.r();
    for(int i=0; i<n; ++i){
        double betx = lattice.betx(i);
        double alfx = lattice.alfx(i);
        double dx = lattice.dx(i);
        double dpx = lattice.dpx(i);
        double alfy = lattice.alfy(i);

        d_tld = alfx*dx+betx*dpx;
        sigma_x = sqrt(sigma_xbet[i]*sigma_xbet[i]+dx*dx*dp_p*dp_p);
        sigma_tmp = dp_p*sigma_xbet[i]/(gamma*sigma_x);
        q = 2*beta*gamma*sqrt(sigma_y[i]/r);

        a_f[i] = sigma_tmp*sqrt(1+alfx*alfx)/sigma_xbetp[i];
        b2_f[i] = sigma_tmp*sqrt(1+alfy*alfy)/sigma_yp[i];
        b2_f[i] *= b2_f[i];
        c2_f[i] = q*sigma_tmp;
        c2_f[i] *= c2_f[i];
        d2_f[i] = dp_p*dx/sigma_x;
        d2_f[i] *= d2_f[i];
        dtld_f[i] = dp_p*d_tld/sigma_x;

        k1_f[i] = 1/c2_f[i];
        k2_f[i] = a_f[i]*a_f[i]*k1_f[i];
        k3_f[i] = b2_f[i]*k1_f[i];
    }
}

void IBSSolver_Martini::coef_f(int nu, int nv)
{
#ifndef NDEBUG
	std::cerr << "ISBSolver_Martini::coef_f()" << std::endl;
#endif
    if(sin_u.get()==nullptr) sin_u.reset(new double[nu]);
    if(sin_u2.get()==nullptr) sin_u2.reset(new double[nu]);
    if(cos_u2.get()==nullptr) cos_u2.reset(new double[nu]);
    if(sin_v.get()==nullptr) sin_v.reset(new double[nv]);
    if(cos_v.get()==nullptr) cos_v.reset(new double[nv]);
    if(sin_u2_cos_v2.get()==nullptr) sin_u2_cos_v2.reset(new double[nu*nv]);

    if(g1.get()==nullptr) g1.reset(new double[nu*nv]);
    if(g2_1.get()==nullptr) g2_1.reset(new double[nu*nv]);
    if(g2_2.get()==nullptr) g2_2.reset(new double[nu*nv]);
    if(g3.get()==nullptr) g3.reset(new double[nu]);

    double du = k_pi/nu;
    double u = -0.5*du;
    for(int i=0; i<nu; ++i){
        u += du;
        sin_u[i] = sin(u);
        sin_u2[i] = sin_u[i]*sin_u[i];
        cos_u2[i] = 1-sin_u2[i];
        g3[i] = 1-3*cos_u2[i];
    }

    double dv = 2*k_pi/nv;
    double v = -0.5*dv;
    for(int i=0; i<nv; ++i){
        v += dv;
        sin_v[i] = sin(v);
        cos_v[i] = cos(v);
    }

    int cnt = 0;
    for(int i=0; i<nu; ++i){
        for(int j=0; j<nv; ++j){
            sin_u2_cos_v2[cnt] = sin_u2[i]*cos_v[j]*cos_v[j];
            g1[cnt] = 1-3*sin_u2_cos_v2[cnt];
            g2_1[cnt] = 1-3*sin_u2[i]*sin_v[j]*sin_v[j];
            g2_2[cnt] = 6*sin_u[i]*sin_v[j]*cos_v[j];
            ++cnt;
        }
    }
}

void IBSSolver_Martini::f(int n_element)
{
    if(f1.get()==nullptr) f1.reset(new double[n_element]);
    if(f2.get()==nullptr) f2.reset(new double[n_element]);
    if(f3.get()==nullptr) f3.reset(new double[n_element]);
    memset(f1.get(), 0, n_element*sizeof(double));
    memset(f2.get(), 0, n_element*sizeof(double));
    memset(f3.get(), 0, n_element*sizeof(double));

if (log_c_ > 0) {
    double duvTimes2Logc = 2*k_pi*k_pi/(nu_*nv_) * 2 * log_c_;
#ifndef NDEBUG
	std::cerr << "ISBSolver_Martini::f(): using log_c" << std::endl;
#endif
//    #pragma omp parallel for
    for(int ie=0; ie<n_element; ++ie){
        const double tmp_a_f = a_f[ie];
        const double tmp_dtld_f = dtld_f[ie];
        const double tmp_b2_f = b2_f[ie];
        const double tmp_k1_f = k1_f[ie];

        for(int iu=0; iu<nu_; ++iu){
            double foo1 = 0, foo2 = 0, foo3 = 0;
            const double tmp_sin_u = sin_u[iu];
            const double tmp_sin_u2 = sin_u2[iu];
            const double tmp_cos_u2 = cos_u2[iu];
            for(int iv=0; iv<nv_; ++iv){
                const int cnt = iu * nv_ + iv;
                double tmp = tmp_a_f * sin_v[iv] - tmp_dtld_f * cos_v[iv];
                tmp *= tmp;
                const double inv_int_z = (sin_u2_cos_v2[cnt] + tmp_sin_u2 * tmp + tmp_b2_f * tmp_cos_u2) * tmp_k1_f;
                foo1 += g1[cnt] / inv_int_z;
                foo2 += (g2_1[cnt] + g2_2[cnt] * tmp_dtld_f / tmp_a_f) / inv_int_z;
                foo3 += g3[iu] / inv_int_z;
            }
            f1[ie] += tmp_sin_u * foo1;
            f2[ie] += tmp_sin_u * foo2;
            f3[ie] += tmp_sin_u * foo3;
        }
        f1[ie] *= k1_f[ie] * duvTimes2Logc;
        f2[ie] *= k2_f[ie] * duvTimes2Logc;
        f3[ie] *= k3_f[ie] * duvTimes2Logc;
    }
} else {
#ifndef NDEBUG
	std::cerr << "ISBSolver_Martini::f(): using nz" << std::endl;
#endif
    for(int ie=0; ie<n_element; ++ie){
        int cnt = 0;
        double duv = 2*k_pi*k_pi/(nv_*nu_);
        for(int iu=0; iu<nu_; ++iu){
            for(int iv=0; iv<nv_; ++iv){
                double tmp = a_f[ie]*sin_v[iv]-dtld_f[ie]*cos_v[iv];
                tmp *= tmp;
                double d_uv = (sin_u2_cos_v2[cnt]+sin_u2[iu]*tmp+b2_f[ie]*cos_u2[iu])*k1_f[ie];
                double int_z = 0;
                double dz = 20/(d_uv*nz_);
                double z = -0.5*dz;
                for(int iz=0; iz<nz_; ++iz){
                    z += dz;
                    int_z += exp(-1*d_uv*z)*log(1+z*z)*dz;
                }
                f1[ie] += sin_u[iu]*int_z*g1[cnt];
                f2[ie] += sin_u[iu]*int_z*(g2_1[cnt]+g2_2[cnt]*dtld_f[ie]/a_f[ie]);
                f3[ie] += sin_u[iu]*int_z*g3[iu];
                ++cnt;
            }
        }
        f1[ie] *= k1_f[ie]*duv;
        f2[ie] *= k2_f[ie]*duv;
        f3[ie] *= k3_f[ie]*duv;
    }
}
}

double IBSSolver_Martini::coef_a(const Lattice &lattice, const Beam &beam) const
{
    double lambda = beam.particle_number()/lattice.circ();
    if(beam.bunched()) lambda = beam.particle_number()/(2*sqrt(k_pi)*beam.sigma_s());

    double beta3 = beam.beta();
    beta3 = beta3*beta3*beta3;
    double gamma4 = beam.gamma();
    gamma4 *= gamma4;
    gamma4 *= gamma4;
    return k_c*beam.r()*beam.r()*lambda/(16*k_pi*sqrt(k_pi)*beam.dp_p()*beta3*gamma4)/(beam.emit_x()*beam.emit_y());
}

void IBSSolver_Martini::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs)
{
    int n_element = lattice.n_element();
    bunch_size(lattice, beam);
    abcdk(lattice, beam);
    if(reset()) {
        coef_f(nu_, nv_);
        reset_off();
    }
    f(n_element);

    double a = coef_a(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    for(int i=0; i<n_element-1; ++i){
        double l_element = lattice.l_element(i);
        rs += n*a*(1-d2_f[i])*f1[i]*l_element;
        rx += a*(f2[i]+f1[i]*(d2_f[i]+dtld_f[i]*dtld_f[i]))*l_element;
        ry += a*f3[i]*l_element;
    }

    double circ = lattice.circ();
    rx /= circ;
    ry /= circ;
    rs /= circ;
    if(k_>0) ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}



IBSSolver_BM::IBSSolver_BM(double log_c, double k)
    : IBSSolver(log_c, k)
{
#ifndef NDEBUG
	std::cerr << "DEBUG: ISBSolver_BM constructor" << std::endl;
#endif
}

void IBSSolver_BM::alloc_var(std::unique_ptr<double []>& ptr, int n) {
    if(ptr.get()==nullptr) {
        ptr = std::unique_ptr<double []>(new double [n]);
    }
}

//Calculate the variables that only depend on the TWISS parameters and energy of the beam.
void IBSSolver_BM::init_fixed_var(const Lattice& lattice, const Beam& beam) {
    int n = lattice.n_element();
    alloc_var(phi, n);
    alloc_var(dx2, n);
    alloc_var(dx_betax_phi_2, n);
    alloc_var(sqrt_betay, n);
    alloc_var(gamma_phi_2, n);

    for(int i=0; i<n; ++i) {
        double dx_betax = lattice.dx(i)/lattice.betx(i);
        dx2[i] = lattice.dx(i) * lattice.dx(i);
        phi[i] = lattice.dpx(i) + lattice.alfx(i)*dx_betax;
        dx_betax_phi_2[i] = dx_betax*dx_betax + phi[i]*phi[i];
        sqrt_betay[i] = sqrt(lattice.bety(i));
        gamma_phi_2[i] = beam.gamma() * beam.gamma() * phi[i] * phi[i];
    }

//    std::ofstream output_particles;
//    output_particles.open("ibs_bm_init_fixed_var.txt");
//    output_particles<<"dx                 "
//                    <<"bet_x              "
//                    <<"dpx                "
//                    <<"alfa_x             "
//                    <<"dx^2               "
//                    <<"phi                "
//                    <<"dx^2/bet_x^2+phi^2 "
//                    <<"sqrt(bet_y)        "
//                    <<"gamma^2*phi^2      "
//        <<std::endl;
//    output_particles.precision(10);
//    output_particles<<std::showpos;
//    output_particles<<std::scientific;
//    for(int i=0; i<n; ++i) {
//        output_particles<<lattice.dx(i)<<' '<<lattice.betx(i)<<' '<<lattice.dpx(i)<<' '<<lattice.alfx(i)<<' '<<dx2[i]<<' '
//        <<phi[i]<<' '<<dx_betax_phi_2[i]<<' '<<sqrt_betay[i]<<' '<<gamma_phi_2[i]<<std::endl;
//    }
//    output_particles.close();


}

void IBSSolver_BM::calc_kernels(const Lattice& lattice, const Beam& beam) {
    auto emit_x = beam.emit_x();
    auto emit_y = beam.emit_y();
    auto sigma_p2 = beam.dp_p() * beam.dp_p();
    auto inv_sigma_p2 = 1/sigma_p2;
    auto n = lattice.n_element();
    auto gamma2 = beam.gamma()*beam.gamma();

//    std::cout<<"emit_x: "<<emit_x<<std::endl
//             <<"emit_y: "<<emit_y<<std::endl
//             <<"gamma:  "<<beam.gamma()<<std::endl;

    alloc_var(psi, n);
    alloc_var(sx, n);
    alloc_var(sp, n);
    alloc_var(sxp, n);
    alloc_var(inv_sigma, n);

//    std::ofstream output_particles;
//    output_particles.open("ibs_bm_calc_kernels.txt");
//    output_particles<<"bet_x              "
//                    <<"bet_y              "
//                    <<"sigma_x            "
//                    <<"sigma_y            "
//                    <<"1/sigma_x/sigma_y  "
//                    <<"a1                 "
//                    <<"a2                 "
//                    <<"lamda_1            "
//                    <<"lamda_2            "
//                    <<"lamda_3            "
//                    <<"psi                "
//                    <<"sxp                "
//                    <<"sp                 "
//                    <<"sx                 "
//        <<std::endl;
//    output_particles.precision(10);
//    output_particles<<std::showpos;
//    output_particles<<std::scientific;





    for(int i=0; i<n; ++i) {
        auto betx = lattice.betx(i);
        auto bety = lattice.bety(i);
        auto sigma_x = sqrt(dx2[i]*sigma_p2+emit_x*betx);
        auto sigma_y = sqrt(emit_y*bety);
//        auto l = lattice.circ();
        inv_sigma[i] = 1/(sigma_x*sigma_y);

        auto ax = betx/emit_x;
        auto lamda_1 = bety/emit_y; //lamda_1 = ay.
        auto as = ax*dx_betax_phi_2[i] + inv_sigma_p2;
        auto a1 = gamma2*as;
        auto a2 = (ax-a1)/2;
        a1 = (ax+a1)/2;

        auto lamda_2 = sqrt(a2*a2+ax*ax*gamma_phi_2[i]);
        auto lamda_3 = a1 - lamda_2;
        auto tmp1 = 3/lamda_2;
        lamda_2 = a1 + lamda_2;

        auto inv_lamda_1 = 1/lamda_1;
        auto inv_lamda_2 = 1/lamda_2;
        auto inv_lamda_3 = 1/lamda_3;

        auto r1 = rd(inv_lamda_2, inv_lamda_3, inv_lamda_1);
        auto r2 = rd(inv_lamda_3, inv_lamda_1, inv_lamda_2);
        auto r3 = 3*sqrt(lamda_1*lamda_2*lamda_3)-r1-r2;

        r1 *= inv_lamda_1*2;
        r2 *= inv_lamda_2;
        r3 *= inv_lamda_3;

        psi[i] = -r1 + r2 + r3;

        sxp[i] = tmp1*ax*gamma_phi_2[i]*(r3-r2);
        tmp1 = tmp1*a2;
        auto tmp2 = 1 + tmp1;
        tmp1 = 1 - tmp1;
        sp[i] = gamma2*(r1-r2*tmp1-r3*tmp2)/2;
        sx[i] = (r1-r2*tmp2-r3*tmp1)/2;


//        output_particles<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<sigma_x<<' '<<sigma_y<<' '<<inv_sigma[i]<<' '
//        <<a1<<' '<<a2<<' '<<lamda_1<<' '<<lamda_2<<' '<<lamda_3<<' '<<psi[i]<<' '<<sxp[i]<<' '<<sp[i]<<' '<<sx[i]<<std::endl;
    }

//    output_particles.close();
}

double IBSSolver_BM::coef_bm(const Lattice &lattice, const Beam &beam) const {
    double lambda = 1;
    if (beam.bunched())
        lambda /= 2*sqrt(k_pi)*beam.sigma_s();
    else
        lambda /= lattice.circ();

    double beta3 = beam.beta()*beam.beta()*beam.beta();
    double gamma5 = beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma();
    return lambda*beam.particle_number()*beam.r()*beam.r()*k_c/(lattice.circ()*6*sqrt(k_pi)*beta3*gamma5);
}

void IBSSolver_BM::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs)
{
    int n_element = lattice.n_element();

    if (reset()) {
        init_fixed_var(lattice, beam);
        reset_off();
    }

    calc_kernels(lattice, beam);
    double c_bm = coef_bm(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    for(int i=0; i<n_element-1; ++i){
        double l_element = lattice.l_element(i);
        rs += inv_sigma[i]*sp[i]*l_element;
        ry += lattice.bety(i)*inv_sigma[i]*psi[i]*l_element;
        rx += lattice.betx(i)*inv_sigma[i]*(sx[i]+dx_betax_phi_2[i]*sp[i]+sxp[i])*l_element;
    }

//    double circ = lattice.circ();
    double lc = log_c();
    c_bm *= lc;
    rs *= n*c_bm;
    rx *= c_bm;
    ry *= c_bm;

    rs /= beam.dp_p()*beam.dp_p();
    rx /= beam.emit_x();
    ry /= beam.emit_y();

    if(k_>0)
        ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}

