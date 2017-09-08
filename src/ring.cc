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
    std::map<std::string, int> keywords = {{"BETX",0}, {"ALFX",0}, {"MUX",0}, {"DX",0}, {"DPX",0}, {"BETY",0}, {"ALFY",0},
                                            {"MUY",0}, {"DY",0}, {"DPY",0}, {"S",0}};


    bool read = false;
    unsigned long int cnt = 0;
    while(std::getline(infile,line)){
        std::istringstream iss(line);
        std::string name;
        std::string keyword;
        double num;
        iss>>name;
//        if (name.compare("\"MACHINE$START\"")==0) read=true;
        std::size_t found = name.find("$START");
        if (found!=std::string::npos) read=true;

        found = name.find("*");
        if (found!=std::string::npos) {
            iss >> name;
            iss >> keyword;
            int cnt_keywords = 0;
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
            assert(cnt_keywords==keywords.size()&&"TWISS PARAMETERS NOT COMPLETE!");
        }

        if(read){
            iss>>keyword;
            int cnt_position = 0;
            int cnt_keywords = 0;
            while(!iss.eof()) {
                iss>>num;
                if(cnt_position == keywords["BETX"]) {betx_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["ALFX"]) {alfx_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["MUX"]) {mux_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["DX"]) {dx_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["DPX"]) {dpx_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["BETY"]) {bety_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["ALFY"]) {alfy_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["MUY"]) {muy_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["DY"]) {dy_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["DPY"]) {dpy_.push_back(num); ++cnt_keywords;}
                else if(cnt_position == keywords["S"]) {s_.push_back(num); ++cnt_keywords;}
                if (cnt_keywords == keywords.size()) {++cnt; break;}
                ++cnt_position;
            }
        }
        found = name.find("$END");
        if (found!=std::string::npos){
//        if(name.compare("\"MACHINE$END\"")==0){
            n_element_ = cnt;
            read = false;
            break;
        }
    }
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
}

Ring::Ring(Lattice &lattice_defined, Beam &beam_defined) {
    beam_ = &beam_defined;
    lattice_ = &lattice_defined;
    circ_ = lattice_->circ();
    if(beam_->bunched())
        beta_s_ = beam_->sigma_s()/beam_->dp_p();
}


