#include <fstream>
#include <memory>
#include "functions.h"
#include "ibs.h"

//Scratch variables for IBS calculation (Martini model)

std::unique_ptr<double []> sigma_xbet, sigma_xbetp, sigma_y, sigma_yp;
std::unique_ptr<double []> a_f, b2_f, c2_f, d2_f, dtld_f, k1_f, k2_f, k3_f;

std::unique_ptr<double []> sin_u, sin_u2, cos_u2, sin_v, cos_v, sin_u2_cos_v2;
std::unique_ptr<double []> g1, g2_1, g2_2, g3;
std::unique_ptr<double []> f1, f2, f3;

//Set ptrs for bunch_size
void set_bunch_size(int n) {
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
int bunch_size(Lattice &lattice, Beam &beam) {
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
    return 0;
}

//Set the ptrs for abcdk
void set_abcdk(int n) {
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
int abcdk(Lattice &lattice, Beam &beam){
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
    return 0;
}

int coef_f(int nu, int nv){
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
    return 0;
}

int f(int n_element, int nu, int nv, int nz){
    if(f1.get()==nullptr) f1.reset(new double[n_element]);
    if(f2.get()==nullptr) f2.reset(new double[n_element]);
    if(f3.get()==nullptr) f3.reset(new double[n_element]);
    memset(f1.get(), 0, n_element*sizeof(double));
    memset(f2.get(), 0, n_element*sizeof(double));
    memset(f3.get(), 0, n_element*sizeof(double));

    for(int ie=0; ie<n_element; ++ie){
        int cnt = 0;
//        f1[ie] = 0;
//        f2[ie] = 0;
//        f3[ie] = 0;
        double duv = 2*k_pi*k_pi/(nv*nu);
        for(int iu=0; iu<nu; ++iu){
            for(int iv=0; iv<nv; ++iv){
                double tmp = a_f[ie]*sin_v[iv]-dtld_f[ie]*cos_v[iv];
                tmp *= tmp;
                double d_uv = (sin_u2_cos_v2[cnt]+sin_u2[iu]*tmp+b2_f[ie]*cos_u2[iu])*k1_f[ie];
                double int_z = 0;
                double dz = 20/(d_uv*nz);
                double z = -0.5*dz;
                for(int iz=0; iz<nz; ++iz){
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
    return 0;
}

int f(int n_element, int nu, int nv, double log_c){
//    if(f1.get()==nullptr) std::cout<<"create f1"<<std::endl;
    if(f1.get()==nullptr) f1.reset(new double[n_element]);
    if(f2.get()==nullptr) f2.reset(new double[n_element]);
    if(f3.get()==nullptr) f3.reset(new double[n_element]);
    memset(f1.get(), 0, n_element*sizeof(double));
    memset(f2.get(), 0, n_element*sizeof(double));
    memset(f3.get(), 0, n_element*sizeof(double));
    for(int ie=0; ie<n_element; ++ie){
        int cnt = 0;
//        f1[ie] = 0;
//        f2[ie] = 0;
//        f3[ie] = 0;
        double duv = 2*k_pi*k_pi/(nu*nv);
        for(int iu=0; iu<nu; ++iu){
            for(int iv=0; iv<nv; ++iv){
                double tmp = a_f[ie]*sin_v[iv]-dtld_f[ie]*cos_v[iv];
                tmp *= tmp;
                double d_uv = (sin_u2_cos_v2[cnt]+sin_u2[iu]*tmp+b2_f[ie]*cos_u2[iu])*k1_f[ie];
                double int_z = 2*log_c/d_uv;
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
    return 0;
}

double coef_a(Lattice &lattice, Beam &beam){
    double lambda = beam.particle_number()/lattice.circ();
    if(beam.bunched()) lambda = beam.particle_number()/(2*sqrt(k_pi)*beam.sigma_s());

    double beta3 = beam.beta();
    beta3 = beta3*beta3*beta3;
    double gamma4 = beam.gamma();
    gamma4 *= gamma4;
    gamma4 *= gamma4;
    return k_c*beam.r()*beam.r()*lambda/(16*k_pi*sqrt(k_pi)*beam.dp_p()*beta3*gamma4)/(beam.emit_x()*beam.emit_y());
}

int ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y){
    double rxc = 0.5*(rx*(2-k)+ry*k*emit_y/emit_x);
    double ryc = 0.5*(ry*(2-k)+rx*k*emit_x/emit_y);
    rx = rxc;
    ry = ryc;
    return 0;
}

int ibs_martini(Lattice &lattice, Beam &beam, IBSParas &ibs_paras,double &rx, double &ry, double &rs) {
    int nu = ibs_paras.nu();
    int nv = ibs_paras.nv();
    int nz = ibs_paras.nz();
    double log_c = ibs_paras.log_c();
    double k = ibs_paras.k();
    int n_element = lattice.n_element();
    bunch_size(lattice, beam);
    abcdk(lattice, beam);
    if(ibs_paras.reset()) {
        coef_f(nu,nv);
        ibs_paras.reset_off();
    }
    if(ibs_paras.use_log_c())
        f(n_element,nu,nv,log_c);
    else
        f(n_element,nu,nv,nz);

    double a = coef_a(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;
    double circ = lattice.circ();

    if(ibs_paras.ibs_by_element()) {
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
            double rsi = n*a*(1-d2_f[i])*f1[i]*l_element*inv_circ;
            double rxi = a*(f2[i]+f1[i]*(d2_f[i]+dtld_f[i]*dtld_f[i]))*l_element*inv_circ;
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
            double l_element = lattice.l_element(i);
            rs += n*a*(1-d2_f[i])*f1[i]*l_element;
            rx += a*(f2[i]+f1[i]*(d2_f[i]+dtld_f[i]*dtld_f[i]))*l_element;
            ry += a*f3[i]*l_element;
        }
        rx /= circ;
        ry /= circ;
        rs /= circ;
    }

    if(k>0) ibs_coupling(rx, ry, k, beam.emit_nx(), beam.emit_ny());
    return 0;
}


// IBS calculation by Bjorken-Mtingwa model using Sergei Nagaitsev's method
std::unique_ptr<double []> phi;
std::unique_ptr<double []> dx2; //D_x * D_x
std::unique_ptr<double []> dx_betax_phi_2; // D_x * D_x / (beta_x * beta_x) + phi * phi
std::unique_ptr<double []> sqrt_betay;  // sqrt(beta_y)
std::unique_ptr<double []> gamma_phi_2; // gamma * gamma * phi * phi
std::unique_ptr<double []> psi;
std::unique_ptr<double []> sx, sp, sxp;
std::unique_ptr<double []> inv_sigma;

void alloc_var(std::unique_ptr<double []>& ptr, int n) {
    if(ptr.get()==nullptr) {
        ptr = std::unique_ptr<double []>(new double [n]);
    }
}

//Calculate the variables that only depend on the TWISS parameters and energy of the beam.
void init_fixed_var(Lattice& lattice, Beam& beam) {
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
}

void calc_kernels(Lattice& lattice, Beam& beam) {
    auto emit_x = beam.emit_x();
    auto emit_y = beam.emit_y();
    auto sigma_p2 = beam.dp_p() * beam.dp_p();
    auto inv_sigma_p2 = 1/sigma_p2;
    auto n = lattice.n_element();
    auto gamma2 = beam.gamma()*beam.gamma();

    alloc_var(psi, n);
    alloc_var(sx, n);
    alloc_var(sp, n);
    alloc_var(sxp, n);
    alloc_var(inv_sigma, n);

    for(int i=0; i<n; ++i) {
        auto betx = lattice.betx(i);
        auto bety = lattice.bety(i);
        auto sigma_x = sqrt(dx2[i]*sigma_p2+emit_x*betx);
        auto sigma_y = sqrt(emit_y*bety);
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

    }
}

double coef_bm(Lattice &lattice, Beam &beam) {
    double lambda = 1;
    if(beam.bunched()) lambda /= 2*sqrt(k_pi)*beam.sigma_s();
    else lambda /= lattice.circ();

    double beta3 = beam.beta()*beam.beta()*beam.beta();
    double gamma5 = beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma();
    return lambda*beam.particle_number()*beam.r()*beam.r()*k_c/(lattice.circ()*6*sqrt(k_pi)*beta3*gamma5);
}

int ibs_bm(Lattice &lattice, Beam &beam, IBSParas &ibs_paras,double &rx, double &ry, double &rs) {
        double k = ibs_paras.k();
    int n_element = lattice.n_element();

    if(ibs_paras.reset()) {
        init_fixed_var(lattice, beam);
        ibs_paras.reset_off();
    }

    calc_kernels(lattice, beam);
    double c_bm = coef_bm(lattice, beam);
    double lc = ibs_paras.log_c();
    c_bm *= lc;

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    if(ibs_paras.ibs_by_element()) {
        std::ofstream out;
        out.open("ibs_by_element.txt");
        out<<"s bet_x bet_y alf_x alf_y disp_x disp_y rx_i ry_i rs_i rx_int ry_int rs_int"<<std::endl;
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            double rxi = (lattice.betx(i)*inv_sigma[i]*(sx[i]+dx_betax_phi_2[i]*sp[i]+sxp[i])*l_element)*c_bm/beam.emit_x();
            double ryi = (lattice.bety(i)*inv_sigma[i]*psi[i]*l_element)*c_bm/beam.emit_y();
            double rsi = (inv_sigma[i]*sp[i]*l_element)*n*c_bm/(beam.dp_p()*beam.dp_p());
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
            rs += inv_sigma[i]*sp[i]*l_element;
            ry += lattice.bety(i)*inv_sigma[i]*psi[i]*l_element;
            rx += lattice.betx(i)*inv_sigma[i]*(sx[i]+dx_betax_phi_2[i]*sp[i]+sxp[i])*l_element;
        }

        rs *= n*c_bm;
        rx *= c_bm;
        ry *= c_bm;

        rs /= beam.dp_p()*beam.dp_p();
        rx /= beam.emit_x();
        ry /= beam.emit_y();
    }

    if(k>0) ibs_coupling(rx, ry, k, beam.emit_nx(), beam.emit_ny());
    return 0;
}

int ibs_rate(Lattice &lattice, Beam &beam, IBSParas &ibs_paras,double &rx, double &ry, double &rs) {
    switch (ibs_paras.model()) {
        case IBSModel::MARTINI : {
            ibs_martini(lattice, beam, ibs_paras, rx, ry, rs);
            break;
        }
        case IBSModel::BM : {
            ibs_bm(lattice, beam, ibs_paras, rx, ry, rs);
            break;
        }
        default : {
            assert(false&&"IBS model NOT exists!");
        }
    }
    return 0;
}
