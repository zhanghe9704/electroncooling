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
    storage_opt.resize(n);

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

        storage_opt[i] = os;
    }
}

void IBSSolver_Martini::coef_f()
{
#ifndef NDEBUG
	std::cerr << "IBSSolver_Martini::coef_f()" << std::endl;
#endif
    storage_u.resize(nu_);
    storage_v.resize(nv_);

    double dv = 2*k_pi/nv_;
    double v = -0.5*dv;
    for(int i=0; i<nv_; ++i){
        v += dv;
        storage_v[i] = TrigonometryStorageV({sin(v), cos(v)});
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
            const double sin_u2_cos_v2 = sin_u2 * storage_v[j].cos_v * storage_v[j].cos_v;
            const double g1 = 1 - 3 * sin_u2_cos_v2;
            const double g2_1 = 1 - 3 * sin_u2 * storage_v[j].sin_v * storage_v[j].sin_v;
            const double g2_2 = 6 * sin_u * storage_v[j].sin_v * storage_v[j].cos_v;
	    TrigonometryStorageUV tempUV = {sin_u2_cos_v2, g1, g2_1, g2_2};
	    uv[j] = tempUV;
        }
        storage_u[i] = TrigonometryStorageU({sin_u, sin_u2, cos_u2, g3, uv});
    }
}

void IBSSolver_Martini::f()
{
    const int n_element = storage_opt.size();
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
            const OpticalStorage &os = storage_opt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu) {
                const TrigonometryStorageU &tu = storage_u[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv) {
                    const TrigonometryStorageV &tv = storage_v[iv];
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
            const OpticalStorage &os = storage_opt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu){
                const TrigonometryStorageU &tu = storage_u[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv){
                    const TrigonometryStorageV &tv = storage_v[iv];
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
    if (cache_invalid) {
        coef_f();
        cache_invalid = false;
    }
    f();

    double a = coef_a(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;
    const double circ = lattice.circ();

    const int n_element = lattice.n_element();
    if(ibs_by_element) {
        std::ofstream out;
        out.open("ibs_by_element.txt");
        out<<"s bet_x bet_y alf_x alf_y disp_x disp_y rx_i ry_i rs_i rx_int ry_int rs_int"<<std::endl;
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        double inv_circ = 1/circ;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            const OpticalStorage &os = storage_opt[i];
            double rsi = n*a*(1-os.d2)*f1[i]*l_element*inv_circ;
            double rxi = a*(f2[i]+f1[i]*(os.d2+os.dtld*os.dtld))*l_element*inv_circ;
            double ryi = a*f3[i]*l_element*inv_circ;
            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();

    }
    else {
        for(int i=0; i<n_element-1; ++i){
            const double l_element = lattice.l_element(i);
            const OpticalStorage &os = storage_opt[i];
            rs += n*a*(1-os.d2)*f1[i]*l_element;
            rx += a*(f2[i]+f1[i]*(os.d2+os.dtld*os.dtld))*l_element;
            ry += a*f3[i]*l_element;
        }
        rx /= circ;
        ry /= circ;
        rs /= circ;
    }

    if(k_>0) ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}



IBSSolver_BM::IBSSolver_BM(double log_c, double k)
    : IBSSolver(log_c, k)
{
#ifndef NDEBUG
	std::cerr << "DEBUG: IBSSolver_BM constructor" << std::endl;
#endif
}

void IBSSolver_BM::init_fixed_var(const Lattice& lattice, const Beam& beam) {
    int n = lattice.n_element();
    optical_strage.resize(n);

    for(int i=0; i<n; ++i) {
        OpticalStorage os;
        double dx_betax = lattice.dx(i)/lattice.betx(i);
        os.dx2 = lattice.dx(i) * lattice.dx(i);
        os.phi = lattice.dpx(i) + lattice.alfx(i)*dx_betax;
        os.dx_betax_phi_2 = dx_betax*dx_betax + os.phi*os.phi;
        os.sqrt_betay = sqrt(lattice.bety(i));
        os.gamma_phi_2 = beam.gamma() * beam.gamma() * os.phi * os.phi;
        optical_strage.at(i)= os;
    }
}

void IBSSolver_BM::calc_kernels(const Lattice& lattice, const Beam& beam) {
    auto emit_x = beam.emit_x();
    auto emit_y = beam.emit_y();
    auto sigma_p2 = beam.dp_p() * beam.dp_p();
    auto inv_sigma_p2 = 1/sigma_p2;
    auto n = lattice.n_element();
    auto gamma2 = beam.gamma()*beam.gamma();

    kernels.resize(n);

    for(int i=0; i<n; ++i) {
        Kernels knl;
        auto betx = lattice.betx(i);
        auto bety = lattice.bety(i);
        auto sigma_x = sqrt(optical_strage.at(i).dx2*sigma_p2+emit_x*betx);
        auto sigma_y = sqrt(emit_y*bety);
        knl.inv_sigma = 1/(sigma_x*sigma_y);

        auto ax = betx/emit_x;
        auto lamda_1 = bety/emit_y; //lamda_1 = ay.
        auto as = ax*optical_strage.at(i).dx_betax_phi_2 + inv_sigma_p2;
        auto a1 = gamma2*as;
        auto a2 = (ax-a1)/2;
        a1 = (ax+a1)/2;

        auto lamda_2 = sqrt(a2*a2+ax*ax*optical_strage.at(i).gamma_phi_2);
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

        knl.psi = -r1 + r2 + r3;

        knl.sxp = tmp1*ax*optical_strage.at(i).gamma_phi_2*(r3-r2);
        tmp1 = tmp1*a2;
        auto tmp2 = 1 + tmp1;
        tmp1 = 1 - tmp1;
        knl.sp = gamma2*(r1-r2*tmp1-r3*tmp2)/2;
        knl.sx = (r1-r2*tmp2-r3*tmp1)/2;
        kernels.at(i) = knl;
    }
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

    if (cache_invalid) {
        init_fixed_var(lattice, beam);
        cache_invalid = false;
    }

    calc_kernels(lattice, beam);
    double c_bm = coef_bm(lattice, beam);
    const double lc = log_c();
    c_bm *= lc;

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    if(ibs_by_element) {

        std::ofstream out;
        out.open("ibs_by_element.txt");
        out<<"s bet_x bet_y alf_x alf_y disp_x disp_y rx_i ry_i rs_i rx_int ry_int rs_int"<<std::endl;
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            double rxi = (lattice.betx(i)*kernels.at(i).inv_sigma*(kernels.at(i).sx
                            +optical_strage.at(i).dx_betax_phi_2*kernels.at(i).sp
                            +kernels.at(i).sxp)*l_element)*c_bm/beam.emit_x();
            double ryi = (lattice.bety(i)*kernels.at(i).inv_sigma*kernels.at(i).psi*l_element)*c_bm/beam.emit_y();
            double rsi = (kernels.at(i).inv_sigma*kernels.at(i).sp*l_element)*n*c_bm/(beam.dp_p()*beam.dp_p());

            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();

    }
    else {
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            rs += kernels.at(i).inv_sigma*kernels.at(i).sp*l_element;
            ry += lattice.bety(i)*kernels.at(i).inv_sigma*kernels.at(i).psi*l_element;
            rx += lattice.betx(i)*kernels.at(i).inv_sigma*(kernels.at(i).sx
                    +optical_strage.at(i).dx_betax_phi_2*kernels.at(i).sp
                    +kernels.at(i).sxp)*l_element;
        }

        rs *= n*c_bm;
        rx *= c_bm;
        ry *= c_bm;

        rs /= beam.dp_p()*beam.dp_p();
        rx /= beam.emit_x();
        ry /= beam.emit_y();
    }

    if(k_>0)
        ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}

