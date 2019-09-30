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
	std::cerr << "DEBUG: IBSSolver_Martini constructor" << std::endl;
#endif
}

//Calculate sigma_xbet, sigma_xbetp, sigma_y, sigma_yp
void IBSSolver_Martini::bunch_size(const Lattice &lattice, const Beam &beam)
{
    int n = lattice.n_element();
    
    sigma_xbet.resize(n);
    sigma_xbetp.resize(n);
    sigma_y.resize(n);
    sigma_yp.resize(n);
    
    double emit_x = beam.emit_x();
    double emit_y = beam.emit_y();
    for(int i=0; i<n; ++i) {
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

//Calculate a, b, c, d and dtld
//Call bunch_size() before calling this one
void IBSSolver_Martini::abcdk(const Lattice &lattice, const Beam &beam)
{
    double d_tld, q, sigma_x, sigma_tmp;
    const int n = lattice.n_element();
    storageOpt.resize(n);

    const double dp_p = beam.dp_p();
    const double beta = beam.beta();
    const double gamma = beam.gamma();
    const double r = beam.r();
    for(int i=0; i<n; ++i){
        const double betx = lattice.betx(i);
        const double alfx = lattice.alfx(i);
        const double dx = lattice.dx(i);
        const double dpx = lattice.dpx(i);
        const double alfy = lattice.alfy(i);

        d_tld = alfx*dx+betx*dpx;
        sigma_x = sqrt(sigma_xbet[i]*sigma_xbet[i]+dx*dx*dp_p*dp_p);
        sigma_tmp = dp_p*sigma_xbet[i]/(gamma*sigma_x);
        q = 2*beta*gamma*sqrt(sigma_y[i]/r);

        OpticalStorage os;
	os.a = sigma_tmp*sqrt(1+alfx*alfx)/sigma_xbetp[i];
        os.b2 = sigma_tmp*sqrt(1+alfy*alfy)/sigma_yp[i];
        os.b2 *= os.b2;
        os.c2 = q*sigma_tmp;
        os.c2 *= os.c2;
        os.d2 = dp_p*dx/sigma_x;
        os.d2 *= os.d2;
        os.dtld = dp_p*d_tld/sigma_x;

        os.k1 = 1.0 / os.c2;
        os.k2 = os.a * os.a * os.k1;
        os.k3 = os.b2 * os.k1;
	
	storageOpt[i] = os;
    }
}

void IBSSolver_Martini::coef_f()
{
#ifndef NDEBUG
	std::cerr << "IBSSolver_Martini::coef_f()" << std::endl;
#endif
    storageU.resize(nu_);
    storageV.resize(nv_);
    
    double dv = 2*k_pi/nv_;
    double v = -0.5*dv;
    for(int i=0; i<nv_; ++i){
        v += dv;
        storageV[i] = TrigonometryStorageV({sin(v), cos(v)});
    }

    double du = k_pi/nu_;
    double u = -0.5*du;
    for(int i=0; i<nu_; ++i){
        u += du;
        double sin_u = sin(u);
        double sin_u2 = sin_u * sin_u;
        double cos_u2 = 1 - sin_u2;
        double g3 = 1 - 3 * cos_u2;
	std::vector<TrigonometryStorageUV> uv;
	uv.resize(nv_);
        for(int j=0; j<nv_; ++j){
            const double sin_u2_cos_v2 = sin_u2 * storageV[j].cos_v * storageV[j].cos_v;
            const double g1 = 1 - 3 * sin_u2_cos_v2;
            const double g2_1 = 1 - 3 * sin_u2 * storageV[j].sin_v * storageV[j].sin_v;
            const double g2_2 = 6 * sin_u * storageV[j].sin_v * storageV[j].cos_v;
	    TrigonometryStorageUV tempUV = {sin_u2_cos_v2, g1, g2_1, g2_2};
	    uv[j] = tempUV;
        }
        storageU[i] = TrigonometryStorageU({sin_u, sin_u2, cos_u2, g3, uv});
    }
}

void IBSSolver_Martini::f()
{
    const int n_element = storageOpt.size();
    f1.resize(n_element);
    f2.resize(n_element);
    f3.resize(n_element);

    if (log_c_ > 0) {
        const double duvTimes2Logc = 2*k_pi*k_pi/(nu_*nv_) * 2 * log_c_;
#ifndef NDEBUG
        std::cerr << "IBSSolver_Martini::f(): using log_c method; nu=" << nu_ << "; nv=" << nv_ << "; log_c=" << log_c_ << std::endl;
#endif
#ifdef _OPENMP
        // Maybe make this runtime-adjustable. SMT gives no benefit here, so choose n = physical cores
        #pragma omp parallel for num_threads(6)
#endif
        for(int ie=0; ie < n_element; ie++) {
            const OpticalStorage &os = storageOpt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu) {
                const TrigonometryStorageU &tu = storageU[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv) {
                    const TrigonometryStorageV &tv = storageV[iv];
                    const TrigonometryStorageUV &tuv = tu.uv[iv];

                    double tmp = os.a * tv.sin_v - os.dtld * tv.cos_v;
                    tmp *= tmp;
                    const double inv_int_z = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os.b2 * tu.cos_u2) * os.k1;
                    sum1 += tuv.g1 / inv_int_z;
                    sum2 += (tuv.g2_1 + tuv.g2_2 * os.dtld / os.a) / inv_int_z;
                    sum3 += tu.g3 / inv_int_z;
                }
                tempf1 += tu.sin_u * sum1;
                tempf2 += tu.sin_u * sum2;
                tempf3 += tu.sin_u * sum3;
            }
            f1[ie] = tempf1 * os.k1 * duvTimes2Logc;
            f2[ie] = tempf2 * os.k2 * duvTimes2Logc;
            f3[ie] = tempf3 * os.k3 * duvTimes2Logc;
        }
    } else {
        const double duv = 2*k_pi*k_pi/(nv_*nu_);
#ifndef NDEBUG
        std::cerr << "IBSSolver_Martini::f(): using nz method; nu=" << nu_ << "; nv=" << nv_ << "; nz=" << nz_ << std::endl;
#endif
#ifdef _OPENMP
        // Maybe make this runtime-adjustable. SMT gives no benefit here, so choose n = physical cores
        #pragma omp parallel for num_threads(6)
#endif
        for(int ie=0; ie<n_element; ++ie){
            const OpticalStorage &os = storageOpt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu){
                const TrigonometryStorageU &tu = storageU[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv){
                    const TrigonometryStorageV &tv = storageV[iv];
                    const TrigonometryStorageUV &tuv = tu.uv[iv];

                    double tmp = os.a * tv.sin_v - os.dtld * tv.cos_v;
                    tmp *= tmp;
                    const double d_uv = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os.b2 * tu.cos_u2) * os.k1;
                    double int_z = 0;
                    const double dz = 20/(d_uv*nz_);
                    double z = -0.5*dz;
                    for(int iz=0; iz<nz_; ++iz){
                        z += dz;
                        int_z += exp(-d_uv*z)*log(1+z*z)*dz;
                    }
                    sum1 += int_z * tuv.g1;
                    sum2 += int_z * (tuv.g2_1 + tuv.g2_2 * os.dtld / os.a);
                    sum3 += int_z * tu.g3;
                }
                tempf1 += tu.sin_u * sum1;
                tempf2 += tu.sin_u * sum2;
                tempf3 += tu.sin_u * sum3;
            }
            f1[ie] = tempf1 * os.k1 * duv;
            f2[ie] = tempf2 * os.k2 * duv;
            f3[ie] = tempf3 * os.k3 * duv;
        }
    }
#ifndef NDEBUG
    std::cerr << "IBSSolver_Martini::f() done" << std::endl;
#endif
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
    bunch_size(lattice, beam);
    abcdk(lattice, beam);
    if (cacheInvalid) {
        coef_f();
        cacheInvalid = false;
    }
    f();

    double a = coef_a(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    const int n_element = lattice.n_element();
    for(int i=0; i<n_element-1; ++i){
        const double l_element = lattice.l_element(i);
        const OpticalStorage &os = storageOpt[i];
        rs += n*a*(1-os.d2)*f1[i]*l_element;
        rx += a*(f2[i]+f1[i]*(os.d2+os.dtld*os.dtld))*l_element;
        ry += a*f3[i]*l_element;
    }

    const double circ = lattice.circ();
    rx /= circ;
    ry /= circ;
    rs /= circ;
    if(k_>0) ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}



IBSSolver_BM::IBSSolver_BM(double log_c, double k)
    : IBSSolver(log_c, k)
{
#ifndef NDEBUG
	std::cerr << "DEBUG: IBSSolver_BM constructor" << std::endl;
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

    if (cacheInvalid) {
        init_fixed_var(lattice, beam);
        cacheInvalid = false;
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

