#include "force.h"
#include <cmath>
#include <cstdio>

#include <iostream>
#include <fstream>

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double *d_perp_e, double *d_paral_e, double time_cooler,
                double *force_tr, double *force_long) {
    double f_const = -4*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
//    double dlt2_eff_e = dlt_paral_e*dlt_paral_e+V2_eff_e;
//    double rho_Lamor = me*1e6*dlt_perp_e/(B*c0*c0);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

    double *dlt2_eff_e = new double[ion_number];
    for(unsigned long int i=0; i<ion_number; ++i){
        dlt2_eff_e[i] = d_paral_e[i]*d_paral_e[i]+v2_eff_e;
    }
    double *rho_lamor = new double[ion_number];
    for(unsigned long int i=0; i<ion_number; ++i) {
        rho_lamor[i] = k_me*1e6*d_perp_e[i]/(magnetic_field*k_c*k_c);
    }
    for(unsigned long int i=0; i<ion_number; ++i){
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = v2+dlt2_eff_e[i];
            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);
            double wp = sqrt(wp_const*density_e[i]);

            //Calculate rho_max
            double rho_max = dlt/wp;
            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
            if(rho_max<rho_max_2) rho_max = rho_max_2;
            double rho_max_3 = dlt*time_cooler;
            if(rho_max>rho_max_3) rho_max = rho_max_3;

            double lc = log((rho_max+rho_min+rho_lamor[i])/(rho_min+rho_lamor[i]));   //Coulomb Logarithm
            //Calculate friction force
            double f = f_const*density_e[i]*lc/(dlt*dlt*dlt);
            force_tr[i] = f*v_tr[i];
            force_long[i] = f*v_long[i];
        }
        else{
            force_tr[i] = 0;
            force_long[i] = 0;
        }
    }
    return 0;
}

int parkhomchuk(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                double temperature, double magnetic_field, double d_perp_e, double d_paral_e, double time_cooler,
                double *force_tr, double *force_long) {
    double f_const = -4*charge_number*charge_number*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = temperature*k_c*k_c/(k_me*1e6);
    double dlt2_eff_e = d_paral_e*d_paral_e+v2_eff_e;
    double rho_lamor = k_me*1e6*d_perp_e/(magnetic_field*k_c*k_c);
    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

    std::cout<<"fixed temp force called!"<<std::endl;

//    std::ofstream out;
//	out.open("force_vtr.txt");
//	out.precision(10);
//    out<<std::showpos;
//    out<<std::scientific;
//    for(int i=0; i<ion_number; ++i)
//        out<<v_tr[i]<<std::endl;
//    out.close();
//    out.open("force_vl.txt");
//	out.precision(10);
//    out<<std::showpos;
//    out<<std::scientific;
//    for(int i=0; i<ion_number; ++i)
//        out<<v_long[i]<<std::endl;
//    out.close();
//    out.open("force_density.txt");
//	out.precision(10);
//    out<<std::showpos;
//    out<<std::scientific;
//    for(int i=0; i<ion_number; ++i)
//        out<<density_e[i]<<std::endl;
//    out.close();

    for(unsigned long int i=0; i<ion_number; ++i){
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = v2+dlt2_eff_e;
            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);
            double wp = sqrt(wp_const*density_e[i]);

            //Calculate rho_max
            double rho_max = dlt/wp;
            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
            if(rho_max<rho_max_2) rho_max = rho_max_2;
            double rho_max_3 = dlt*time_cooler;
            if(rho_max>rho_max_3) rho_max = rho_max_3;

            double lc = log((rho_max+rho_min+rho_lamor)/(rho_min+rho_lamor));   //Coulomb Logarithm
            //Calculate friction force
            double f = f_const*density_e[i]*lc/(dlt*dlt*dlt);
            force_tr[i] = f*v_tr[i];
            force_long[i] = f*v_long[i];
        }
        else{
            force_tr[i] = 0;
            force_long[i] = 0;
        }
    }
    return 0;
}

int friction_force(int charge_number, unsigned long int ion_number, double *v_tr, double *v_long, double *density_e,
                   ForceParas &force_paras, double *force_tr, double *force_long){
    switch (force_paras.formula()) {
        case ForceFormula::PARKHOMCHUK: {
            double temperature = force_paras.park_temperature_eff();
            double time_cooler = force_paras.time_cooler();
            double magnetic_field = force_paras.magnetic_field();
            if(force_paras.ptr_d_perp_e()||force_paras.ptr_d_paral_e()) {
                double *d_perp_e = force_paras.ptr_d_perp_e();
                double *d_paral_e = force_paras.ptr_d_paral_e();
                parkhomchuk(charge_number, ion_number, v_tr, v_long, density_e, temperature, magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long);
            }
            else {
                double d_perp_e = force_paras.d_perp_e();
                double d_paral_e = force_paras.d_paral_e();
                parkhomchuk(charge_number, ion_number, v_tr, v_long, density_e, temperature,  magnetic_field, d_perp_e,
                            d_paral_e, time_cooler, force_tr, force_long);
            }
            break;
        }
        default: {
            perror("Choose your formula for friction force calculation!");
            break;
        }
    }
    return 0;
}

