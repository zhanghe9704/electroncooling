#include "ring.h"

#include <assert.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>


Lattice::Lattice(std::string filename) {
    std::ifstream infile;
    infile.open(filename.c_str());
    if(!infile) std::cout<<"Error: failed to load the lattice file!"<<std::endl;
    std::string line;
    std::map<std::string, int> keywords = {{"BETX",10}, {"ALFX",10}, {"DX",10}, {"DPX",10}, {"BETY",10}, {"ALFY",10},
                                           {"DY",10}, {"DPY",10}, {"S",10}};

    std::map<int, std::string> sorted_key;
    bool read = false;
    unsigned long int cnt = 0;
    int start_reading = -2;
    while(std::getline(infile,line)){
        if (line.empty()) continue;
        if (infile.eof()) break;
        std::istringstream iss(line);
        std::string name;
        std::string keyword;
        double num;
        if (!read) ++start_reading;
        if (2==start_reading) read = true;  //Start reading data from the second line after "*..." .
        if (start_reading<1) start_reading = -1;
        iss>>name;

        std::size_t found = name.find("*");
        if (found!=std::string::npos) {
            start_reading = 0;
            iss >> name;
//            iss >> keyword;
            unsigned int cnt_keywords = 0;
            int cnt_position = 0;
            while(cnt_keywords!=keywords.size()) {
                iss >> keyword;
                if (iss.eof()) {
                    break;
                }
                else {
                    auto iter = keywords.find(keyword);
                    if (iter != keywords.end()) {
                        iter->second = cnt_position;
                        ++cnt_keywords;
                    }
                }
                ++cnt_position;
            }
            if (cnt_keywords==keywords.size()){} //Do nothing.
            else if(cnt_keywords==keywords.size()-2 && keywords["DY"]==10 &&
                    keywords["DPY"]==10){} //Do nothing;
            else
                assert(false&&"TWISS PARAMETERS NOT COMPLETE!");

            for (auto iter=keywords.begin(); iter!=keywords.end(); ++iter)
                sorted_key.insert({iter->second, iter->first});
        }

        if(read){
            double s, betx, bety, alfx,alfy, dx, dy, dpx, dpy;
            int cnt_position = 0;
            for (auto iter=sorted_key.begin(); iter!=sorted_key.end(); ++iter) {
                while (cnt_position<iter->first) {
                    iss>>keyword;
                    ++cnt_position;
                }
                iss>>keyword;
                try {
                  num = std::stod(keyword.c_str());
                }
                catch (...) {
                  num = 0;
                }
                if(cnt_position == keywords["BETX"]) {betx = num;}
                else if(cnt_position == keywords["ALFX"]) {alfx = num;}
                else if(cnt_position == keywords["DX"]) {dx = num;}
                else if(cnt_position == keywords["DPX"]) {dpx = num;}
                else if(cnt_position == keywords["BETY"]) {bety = num;}
                else if(cnt_position == keywords["ALFY"]) {alfy = num;}
                else if(cnt_position == keywords["DY"]) {dy = num;}
                else if(cnt_position == keywords["DPY"]) {dpy = num;}
                else if(cnt_position == keywords["S"]) {s = num;}
                else {assert(false&&"Error in lattice file parsing!");}
                ++cnt_position;
            }
            if (cnt>0 && s_.back()==s) continue;

            if (keywords["DY"]==10) dy = 0;
            if (keywords["DPY"]==10) dpy = 0;
            betx_.push_back(betx);
            alfx_.push_back(alfx);
            dx_.push_back(dx);
            dpx_.push_back(dpx);
            bety_.push_back(bety);
            alfy_.push_back(alfy);
            dy_.push_back(dy);
            dpy_.push_back(dpy);
            s_.push_back(s);
            ++cnt;
        }
    }
    n_element_ = cnt;
    circ_ = s_.at(n_element_-1) - s_.at(0);
    for(int i=0; i<n_element_-1; ++i){
        l_element_.push_back(s_.at(i+1)-s_.at(i));
    }
    l_element_.push_back(l_element_.at(0));
}

Ring::Ring(double circ, Beam &beam_defined):circ_(circ) {
    beam_ = &beam_defined;
    lattice_ = nullptr;
    if(beam_->bunched())
        beta_s_ = beam_->sigma_s()/beam_->dp_p();
    f0_ = beam_->beta()*k_c/circ_;
    w0_ = 2*k_pi*f0_;
//    if(!rf) {
//        if(rf->gamma_tr>0) {
//            double gamma = beam_->gamma();
//            slip_factor_ = 1/(rf->gamma_tr*rf->gamma_tr) - 1/(gamma*gamma);
//        }
//    }
}

Ring::Ring(Lattice &lattice_defined, Beam &beam_defined) {
    beam_ = &beam_defined;
    lattice_ = &lattice_defined;
    circ_ = lattice_->circ();
    if(beam_->bunched())
        beta_s_ = beam_->sigma_s()/beam_->dp_p();
    f0_ = beam_->beta()*k_c/circ_;
    w0_ = 2*k_pi*f0_;
//    if(!rf) {
//        if(rf->gamma_tr>0) {
//            double gamma = beam_->gamma();
//            slip_factor_ = 1/(rf->gamma_tr*rf->gamma_tr) - 1/(gamma*gamma);
//        }
//    }
}

void Ring::set_rf() {
    if(rf.gamma_tr>0) {
        double gamma = beam_->gamma();
        slip_factor_ = 1/(rf.gamma_tr*rf.gamma_tr) - 1/(gamma*gamma);
    }
}

double Ring::calc_sync_tune_by_rf() {
    assert(rf.v>0&&rf.gamma_tr>0&&"DEFINE THE RF CAVITY FOR SYNCHROTRON TUNE CALCULATION!"); //RF has to be defined.

    double cp = beam_->gamma()*beam_->beta()*beam_->mass()*1e6; //cp in [eV].
    double tune = rf.h*slip_factor_*rf.v*cos(rf.phi)/(2*k_pi*beam_->beta()*cp);
    if(tune<0) {
        tune = sqrt(-tune);
    }
    else {
        tune = sqrt(tune);
    }

    return tune;
}

double Ring::calc_rf_voltage() {
    if(beam_->bunched()) {
        double energy = beam_->gamma()*beam_->mass()*1e6; //Total energy in [eV].
        double v_rf = 2*k_pi*k_c*k_c*beam_->beta()*beam_->beta()*beam_->beta()*beam_->beta()*slip_factor_*energy*
            beam_->dp_p()*beam_->dp_p()/(w0_*w0_*rf.h*cos(rf.phi)*beam_->sigma_s()*beam_->sigma_s());
        return (v_rf>0)?v_rf:-v_rf;
    }
    else {
        return 0;
    }
}


