#include "ring.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

Lattice::Lattice(std::string filename) {
    std::ifstream infile;
    infile.open(filename.c_str());
    if(!infile) std::cout<<"Error: failed to load the lattice file!"<<std::endl;
    std::string line;

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

        if(read){
            iss>>keyword;
            iss>>num;
            if(s_.empty()||num>s_.at(cnt-1)){
                s_.push_back(num);
                iss>>num;
                betx_.push_back(num);
                iss>>num;
                alfx_.push_back(num);
                iss>>num;
                mux_.push_back(num);
                iss>>num;
                dx_.push_back(num);
                iss>>num;
                dpx_.push_back(num);
                iss>>num;
                bety_.push_back(num);
                iss>>num;
                alfy_.push_back(num);
                iss>>num;
                muy_.push_back(num);
                iss>>num;
                dy_.push_back(num);
                iss>>num;
                dpy_.push_back(num);
                ++cnt;
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


