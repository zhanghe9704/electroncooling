
#include "arbitrary_electron_beam.h"
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

//Load electrons from a given file

long int read_binary_file(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<double> &vx,
                     std::vector<double> &vy, std::vector<double> &vz, std::string filename, long int n = 0, long int skip = 0,
                     int n_buffer = 1000);

long int read_ascii_file(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
                     std::vector<double>& vy, std::vector<double>& vz, std::string filename, long int n = 0,
                     int line_skip = 0);

unsigned long int load_electrons(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
                                 std::vector<double>& vy, std::vector<double>& vz, std::string filename, long int n,
                                int skip, bool binary, int n_buffer) {
    long int n_loaded = 0;
    if (binary) {
        n_loaded = read_binary_file(x, y, z, vx, vy, vz, filename, n, skip, n_buffer);
    }
    else {
        n_loaded = read_ascii_file(x, y, z, vx, vy, vz, filename, n, skip);
    }

    return n_loaded;

}

////Load electrons from a given file
//unsigned long int load_electrons(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
//                                 std::vector<double>& vy, std::vector<double>& vz, unsigned long int n, std::string filename,
//                                int line_skip){
//    std::ifstream infile;
//    infile.open(filename.c_str());
//    if(!infile) assert(false&&"Error: fail to load the electrons from the given file!");
//
//    std::string line;
//    for(int i=0; i<line_skip; ++i) std::getline(infile, line);
//
//    unsigned long int n_loaded = 0;
//    std::string val_str;
//    double val;
//    int j = 0;
//    while(n_loaded<n && std::getline(infile,line)) {
//        if (line.empty()) continue;
//        if (infile.eof()) break;
//        std::istringstream iss(line);
//        for(int i=0; i<6; ++i) {
//            iss>>val_str;
//            try {
//              val = std::stod(val_str.c_str());
//            }
//            catch (...) {
//              val = 0;
//            }
//            switch (i) {
//            case 0: {
//                x.push_back(val);
//                break;
//            }
//            case 1: {
//                y.push_back(val);
//                break;
//            }
//            case 2: {
//                z.push_back(val);
//                break;
//            }
//            case 3: {
//                vx.push_back(val);
//                break;
//            }
//            case 4: {
//                vy.push_back(val);
//                break;
//            }
//            case 5: {
//                vz.push_back(val);
//                break;
//            }
//            }
//        }
//        ++n_loaded;
//        ++j;
//    }
//    return n_loaded;
//}


long int read_ascii_file(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
                                 std::vector<double>& vy, std::vector<double>& vz, std::string filename, long int n,
                                int line_skip){
    std::ifstream infile;
    infile.open(filename.c_str());
    if(!infile) assert(false&&"Error in reading 6D coordinates from ascii file: fail to open the given file!");

    std::string line;
    std::getline(infile,line);
    std::istringstream iss(line);
    long int n_line = 0;
    std::string val_str;
    iss>>val_str;
    try {
        n_line = std::stol(val_str.c_str());
    }
    catch (...) {
        assert(false&&"Error in reading 6D coordinates from ascii file: the 1st line should be the number of particle coordinates stored in the file (an integer) !");
    }


    if (n>0) {
        assert(n<=n_line&&"Error in reading electron 6D coordinates from ascii file: No enough data!");
        n_line = n;
    }

    x.clear();
    y.clear();
    z.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    x.reserve(n_line);
    y.reserve(n_line);
    z.reserve(n_line);
    vx.reserve(n_line);
    vy.reserve(n_line);
    vz.reserve(n_line);

    for(int i=0; i<line_skip; ++i) std::getline(infile, line);

    long int n_loaded = 0;
    double val;
    int j = 0;
    while(n_loaded<n_line && std::getline(infile,line)) {
        if (line.empty()) continue;
        if (infile.eof()) break;
        std::istringstream iss(line);
        for(int i=0; i<6; ++i) {
            iss>>val_str;
            try {
              val = std::stod(val_str.c_str());
            }
            catch (...) {
              val = 0;
            }
            switch (i) {
            case 0: {
                x.push_back(val);
                break;
            }
            case 1: {
                y.push_back(val);
                break;
            }
            case 2: {
                z.push_back(val);
                break;
            }
            case 3: {
                vx.push_back(val);
                break;
            }
            case 4: {
                vy.push_back(val);
                break;
            }
            case 5: {
                vz.push_back(val);
                break;
            }
            }
        }
        ++n_loaded;
        ++j;
    }
    return n_loaded;
}

void read_buffer(std::ifstream &binary_file, std::vector<double> &buffer, std::vector<double> &x, std::vector<double> &y,
                 std::vector<double> &z, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz) {
    binary_file.read((char*)buffer.data(), buffer.size()*sizeof(double));
    assert(binary_file&&"Error in reading electron 6D coordinates from binary file!");
    for(auto itr=buffer.begin(); itr!=buffer.end();) {
        x.push_back(*itr);
        itr++;
        y.push_back(*itr);
        itr++;
        z.push_back(*itr);
        itr++;
        vx.push_back(*itr);
        itr++;
        vy.push_back(*itr);
        itr++;
        vz.push_back(*itr);
        itr++;
    }
}


long int read_binary_file(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<double> &vx,
                     std::vector<double> &vy, std::vector<double> &vz, std::string filename, long int n, long int skip, int n_buffer) {
    std::ifstream binary_file(filename.c_str(), std::ios::in | std::ios::binary);
    long int n_row = 0;
    binary_file.read((char*)&n_row,4);  //Read the number of rows, signed long it: 32 bits/4 bytes
    assert(binary_file&&"Error in reading electron 6D coordinates from binary file: reading failed!");

    binary_file.ignore(skip*6*sizeof(double));
    assert(n+skip<=n_row&&"Error in reading electron 6D coordinates from binary file: No enough data!");
    if (n>0) n_row = n;

    x.clear();
    y.clear();
    z.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    x.reserve(n_row);
    y.reserve(n_row);
    z.reserve(n_row);
    vx.reserve(n_row);
    vy.reserve(n_row);
    vz.reserve(n_row);

    int n_loop = n_row/n_buffer;
    int n_last_chunk = n_row%n_buffer;
    int n_column = 6;

    std::vector<double> buffer(n_buffer*n_column, 0);

    for(int i=0; i<n_loop; ++i)
        read_buffer(binary_file, buffer, x, y, z, vx, vy, vz);
    if(n_last_chunk>0) {
        buffer.resize(n_last_chunk*n_column);
        read_buffer(binary_file, buffer, x, y, z, vx, vy, vz);
    }

    return n_row;

}



//Find the center and the box size of the root box
int find_root_center(double *x, double * y, double * z, const unsigned long int N, double &cx, double &cy, double &cz,
                     double &size) {
	double max_x = x[0];
	double max_y = y[0];
	double max_z = z[0];
	double min_x = x[0];
	double min_y = y[0];
	double min_z = z[0];
	for(unsigned long int i=1; i<N; ++i){
		if (max_x<x[i]) max_x = x[i];
		if (max_y<y[i]) max_y = y[i];
		if (max_z<z[i]) max_z = z[i];
		if (min_x>x[i]) min_x = x[i];
		if (min_y>y[i]) min_y = y[i];
		if (min_z>z[i]) min_z = z[i];
	}
	cx = 0.5*(max_x+min_x);
	cy = 0.5*(max_y+min_y);
	cz = 0.5*(max_z+min_z);
	size = max_x-min_x;
	if (size<(max_y-min_y)) size = max_y - min_y;
	if (size<(max_z-min_z)) size = max_z - min_z;
	size *= (1+1e-6);

	return 0;
}

int find_root_center(std::vector<double>& x, std::vector<double>&  y, std::vector<double>&  z, const unsigned long int n,
                     double &cx, double &cy, double &cz, double &size) {
	double max_x = x[0];
	double max_y = y[0];
	double max_z = z[0];
	double min_x = x[0];
	double min_y = y[0];
	double min_z = z[0];
	for(unsigned long int i=1; i<n; ++i){
		if (max_x<x[i]) max_x = x[i];
		if (max_y<y[i]) max_y = y[i];
		if (max_z<z[i]) max_z = z[i];
		if (min_x>x[i]) min_x = x[i];
		if (min_y>y[i]) min_y = y[i];
		if (min_z>z[i]) min_z = z[i];
	}
	cx = 0.5*(max_x+min_x);
	cy = 0.5*(max_y+min_y);
	cz = 0.5*(max_z+min_z);
	size = max_x-min_x;
	if (size<(max_y-min_y)) size = max_y - min_y;
	if (size<(max_z-min_z)) size = max_z - min_z;
	size *= (1+1e-6);

	return 0;
}


//Create the hierarchical tree structure of the boxes
int create_e_tree(double * x, double * y, double * z, const unsigned long int n, const unsigned int s, vector<Box> &tree,
                unsigned long int * list) {

	Box empty_box;
	double cx, cy, cz, size;	//the center and the size of a box

	//create the root box
	tree.push_back(empty_box);
	find_root_center(x, y, z, n, cx, cy, cz, size);
	tree[0].center[0] = cx;
	tree[0].center[1] = cy;
	tree[0].center[2] = cz;
	tree[0].box_size = size;
	tree[0].n_ptcl = n;

	//prepare for the particle list
	for(unsigned long int i=0; i<n; ++i) {
		list[i] = i+1;
	}
	unsigned long int tmp_idx[8] = {n, n, n, n, n, n, n, n};
	for (unsigned long int i=0; i!=tree.size();++i) {
		if (tree[i].n_ptcl>s) {
			//split the box
			//go through all the particles in the box
			unsigned long int idx=tree[i].first_ptcl;
			while (idx!=n) {	//when idx==n, reached the end of the particle list
				int nx = (x[idx]-tree[i].center[0]>0)?1:0;
				int ny = (y[idx]-tree[i].center[1]>0)?1:0;
				int nz = (z[idx]-tree[i].center[2]>0)?1:0;

				int idx_cld = 4*nx+2*ny+nz;
				if (tree[i].child[idx_cld]>0) {	//child box has created.
					//update the child box
					unsigned long int pos_cld = tree[i].child[idx_cld];
					tree[pos_cld].n_ptcl += 1;
					//udpate the linked list of the particles
					list[tmp_idx[idx_cld]] = idx;
					tmp_idx[idx_cld] = idx;
				}
				else{	//child box has NOT created.
					//create the child box
					unsigned long int pos_cld = tree.size();
					tree.push_back(empty_box);
					tree[pos_cld].parent = i;
					tree[pos_cld].box_size = 0.5*tree[i].box_size;
					tree[pos_cld].first_ptcl = idx;
					tree[pos_cld].n_ptcl +=1;

					if (nx==0) nx=-1;
					if (ny==0) ny=-1;
					if (nz==0) nz=-1;
					tree[pos_cld].center[0] = tree[i].center[0]+0.5*nx*tree[pos_cld].box_size;
					tree[pos_cld].center[1] = tree[i].center[1]+0.5*ny*tree[pos_cld].box_size;
					tree[pos_cld].center[2] = tree[i].center[2]+0.5*nz*tree[pos_cld].box_size;

					tmp_idx[idx_cld] = idx;	//record the last particle in the box

					//update the parent box
					tree[i].child[idx_cld] = pos_cld;
					tree[i].n_child += 1;
					tree[i].first_ptcl = n;
				}
				idx = list[idx];
			}
			tree[i].n_ptcl = 0;
			for (int i=0; i<8; ++i) {
				if (tmp_idx[i]<n){
					list[tmp_idx[i]] = n;
					tmp_idx[i] = n;
				}
			}
		}
	}
//	//Make sure the index of the child boxes stored in the first n_child elements of the child array.
//	for(unsigned long int i=0; i!=tree.size(); ++i) {
//		if (tree[i].n_child>0) {
//			int k=7;
//			for(int j=0; j<tree[i].n_child; ++j) {
//				if (tree[i].child[j]==0) {
//					for(;k>j;--k){
//						if (tree[i].child[k]>0) {
//							std::swap(tree[i].child[j], tree[i].child[k]);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
	return 0;
}

//Create the hierarchical tree structure of the boxes
int create_e_tree(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const unsigned long int n,
                  const unsigned int s, vector<Box> &tree, std::vector<unsigned long int>& list) {

	Box empty_box;
	double cx, cy, cz, size;	//the center and the size of a box

	//create the root box
	tree.push_back(empty_box);
	find_root_center(x, y, z, n, cx, cy, cz, size);
	tree[0].center[0] = cx;
	tree[0].center[1] = cy;
	tree[0].center[2] = cz;
	tree[0].box_size = size;
	tree[0].n_ptcl = n;

	//prepare for the particle list
	list.clear();
	list.reserve(n);
	for(unsigned long int i=0; i<n; ++i) {
        list.push_back(i+1);
//		list[i] = i+1;
	}
	unsigned long int tmp_idx[8] = {n, n, n, n, n, n, n, n};
	for (unsigned long int i=0; i!=tree.size();++i) {
		if (tree[i].n_ptcl>s) {
			//split the box
			//go through all the particles in the box
			unsigned long int idx=tree[i].first_ptcl;
			while (idx!=n) {	//when idx==n, reached the end of the particle list
				int nx = (x.at(idx)-tree[i].center[0]>0)?1:0;
				int ny = (y.at(idx)-tree[i].center[1]>0)?1:0;
				int nz = (z.at(idx)-tree[i].center[2]>0)?1:0;

				int idx_cld = 4*nx+2*ny+nz;
				if (tree[i].child[idx_cld]>0) {	//child box has created.
					//update the child box
					unsigned long int pos_cld = tree[i].child[idx_cld];
					tree[pos_cld].n_ptcl += 1;
					//udpate the linked list of the particles
					list.at(tmp_idx[idx_cld]) = idx;
					tmp_idx[idx_cld] = idx;
				}
				else{	//child box has NOT created.
					//create the child box
					unsigned long int pos_cld = tree.size();
					tree.push_back(empty_box);
					tree[pos_cld].parent = i;
					tree[pos_cld].box_size = 0.5*tree[i].box_size;
					tree[pos_cld].first_ptcl = idx;
					tree[pos_cld].n_ptcl +=1;

					if (nx==0) nx=-1;
					if (ny==0) ny=-1;
					if (nz==0) nz=-1;
					tree[pos_cld].center[0] = tree[i].center[0]+0.5*nx*tree[pos_cld].box_size;
					tree[pos_cld].center[1] = tree[i].center[1]+0.5*ny*tree[pos_cld].box_size;
					tree[pos_cld].center[2] = tree[i].center[2]+0.5*nz*tree[pos_cld].box_size;

					tmp_idx[idx_cld] = idx;	//record the last particle in the box

					//update the parent box
					tree[i].child[idx_cld] = pos_cld;
					tree[i].n_child += 1;
					tree[i].first_ptcl = n;
				}
				idx = list.at(idx);
			}
			tree[i].n_ptcl = 0;
			for (int i=0; i<8; ++i) {
				if (tmp_idx[i]<n){
					list.at(tmp_idx[i]) = n;
					tmp_idx[i] = n;
				}
			}
		}
	}
	return 0;
}



//Create the hierarchical tree structure of the boxes
int create_ion_tree(double * x, double * y, double * z, const unsigned long int n, vector<Box> &tree,
                unsigned long int * list, unsigned long int &idx_out) {

//	Box empty_box;
//	double cx, cy, cz, size;	//the center and the size of a box
//
//	//create the root box
//	tree.push_back(empty_box);
//	find_root_center(x, y, z, n, cx, cy, cz, size);
//	tree[0].center[0] = cx;
//	tree[0].center[1] = cy;
//	tree[0].center[2] = cz;
//	tree[0].box_size = size;
//	tree[0].n_ptcl = n;

	//prepare for the particle list
	unsigned long int idx_in = n;
	idx_out= n;
	double cx, cy, cz, half_box;
    cx = tree[0].center[0];
    cy = tree[0].center[1];
    cz = tree[0].center[2];
    half_box = 0.5*tree[0].box_size;
    unsigned long int ptr_in = n;
    unsigned long int ptr_out = n;
    unsigned long int cnt_in = 0;
	for(unsigned long int i=0; i<n; ++i) {
        bool x_in = fabs(x[i]-cx)<=half_box;
        bool y_in = fabs(y[i]-cy)<=half_box;
        bool z_in = fabs(z[i]-cz)<=half_box;
        if(x_in&&y_in&&z_in) { //ion particle is inside the box
            ++cnt_in;
            if(idx_in == n) idx_in = i;
            if(ptr_in == n) {
                ptr_in = i;
                list[i] = n;
            }
            else {
                list[ptr_in] = i;
                ptr_in = i;
                list[i] = n;
            }
        }
        else {                  //ion particle is outside the box
            if(idx_out == n) idx_out = i;
            if(ptr_out == n) {
                ptr_out = i;
                list[i] = n;
            }
            else {
                list[ptr_out] = i;
                ptr_out = i;
                list[i] = n;
            }
        }
	}
	tree[0].first_ion = idx_in;
	tree[0].n_ion = cnt_in;
    unsigned long int tmp_idx[8] = {n, n, n, n, n, n, n, n};
	for (unsigned long int i=0; i!=tree.size();++i) {
        if(tree[i].n_child>0) {
            for(auto cld_id : tree[i].child) {
                if (cld_id>0) tree[cld_id].first_ion = n;
            }
            unsigned long int idx=tree[i].first_ion;
            unsigned int n_ion = tree[i].n_ion;
            while (n_ion>0&&idx!=n) { //when idx==n, reached the end of the ion list
                int nx = (x[idx]-tree[i].center[0]>0)?1:0;
				int ny = (y[idx]-tree[i].center[1]>0)?1:0;
				int nz = (z[idx]-tree[i].center[2]>0)?1:0;

				int idx_cld = 4*nx+2*ny+nz;
				if (tree[i].child[idx_cld]>0) {	//child box has created --> ion inside electron bunch.
					//update the child box
					unsigned long int pos_cld = tree[i].child[idx_cld];
					tree[pos_cld].n_ion += 1;
					//udpate the linked list of the particles
					if(tmp_idx[idx_cld] == n) {
                        tmp_idx[idx_cld] = idx;
                        tree[pos_cld].first_ion = idx;
                        tree[i].first_ion = n;
					}
					else {
                        list[tmp_idx[idx_cld]] = idx;
                        tmp_idx[idx_cld] = idx;
					}
				}
				else{	//child box has NOT created --> ion is outside electron bunch.
				    if(idx_out == n) idx_out = idx;
				    if(ptr_out == n) {
                        ptr_out = idx;
                    }
                    else {
                        list[ptr_out] = idx;
                        ptr_out = idx;
                    }
				}
				idx = list[idx];
            }
            tree[i].n_ion = 0;
			for (int i=0; i<8; ++i) {
				if (tmp_idx[i]<n){
					list[tmp_idx[i]] = n;
					tmp_idx[i] = n;
				}
			}
        }
	}
	list[ptr_out] = n;

//	//Make sure the index of the child boxes stored in the first n_child elements of the child array.
//	for(unsigned long int i=0; i!=tree.size(); ++i) {
//		if (tree[i].n_child>0) {
//			int k=7;
//			for(int j=0; j<tree[i].n_child; ++j) {
//				if (tree[i].child[j]==0) {
//					for(;k>j;--k){
//						if (tree[i].child[k]>0) {
//							std::swap(tree[i].child[j], tree[i].child[k]);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
	return 0;
}

int create_ion_tree(double * x, double * y, double * z, const unsigned int n, vector<Box> &tree,
                std::vector<unsigned int>& list, unsigned int &idx_out) {
	//prepare for the particle list
	unsigned long int idx_in = n;
	idx_out= n;
	double cx, cy, cz, half_box;
    cx = tree[0].center[0];
    cy = tree[0].center[1];
    cz = tree[0].center[2];
    half_box = 0.5*tree[0].box_size;
    unsigned long int ptr_in = n;
    unsigned long int ptr_out = n;
    unsigned long int cnt_in = 0;
    list.resize(n, 0);
	for(unsigned long int i=0; i<n; ++i) {
        bool x_in = fabs(x[i]-cx)<=half_box;
        bool y_in = fabs(y[i]-cy)<=half_box;
        bool z_in = fabs(z[i]-cz)<=half_box;
        if(x_in&&y_in&&z_in) { //ion particle is inside the box
            ++cnt_in;
            if(idx_in == n) idx_in = i;
            if(ptr_in == n) {
                ptr_in = i;
                list.at(i) = n;
            }
            else {
                list.at(ptr_in) = i;
                ptr_in = i;
                list.at(i) = n;
            }
        }
        else {                  //ion particle is outside the box
            if(idx_out == n) idx_out = i;
            if(ptr_out == n) {
                ptr_out = i;
                list.at(i) = n;
            }
            else {
                list.at(ptr_out) = i;
                ptr_out = i;
                list.at(i) = n;
            }
        }
	}
	tree[0].first_ion = idx_in;
	tree[0].n_ion = cnt_in;
    unsigned long int tmp_idx[8] = {n, n, n, n, n, n, n, n};
	for (unsigned long int i=0; i!=tree.size();++i) {
        if(tree[i].n_child>0) {
            for(auto cld_id : tree[i].child) {
                if (cld_id>0) tree[cld_id].first_ion = n;
            }
            unsigned long int idx=tree[i].first_ion;
            unsigned int n_ion = tree[i].n_ion;
            while (n_ion>0&&idx!=n) { //when idx==n, reached the end of the ion list
                int nx = (x[idx]-tree[i].center[0]>0)?1:0;
				int ny = (y[idx]-tree[i].center[1]>0)?1:0;
				int nz = (z[idx]-tree[i].center[2]>0)?1:0;

				int idx_cld = 4*nx+2*ny+nz;
				if (tree[i].child[idx_cld]>0) {	//child box has created --> ion inside electron bunch.
					//update the child box
					unsigned long int pos_cld = tree[i].child[idx_cld];
					tree[pos_cld].n_ion += 1;
					//udpate the linked list of the particles
					if(tmp_idx[idx_cld] == n) {
                        tmp_idx[idx_cld] = idx;
                        tree[pos_cld].first_ion = idx;
                        tree[i].first_ion = n;
					}
					else {
                        list.at(tmp_idx[idx_cld]) = idx;
                        tmp_idx[idx_cld] = idx;
					}
				}
				else{	//child box has NOT created --> ion is outside electron bunch.
				    if(idx_out == n) idx_out = idx;
				    if(ptr_out == n) {
                        ptr_out = idx;
                    }
                    else {
                        list.at(ptr_out) = idx;
                        ptr_out = idx;
                    }
				}
				idx = list[idx];
            }
            tree[i].n_ion = 0;
			for (int i=0; i<8; ++i) {
				if (tmp_idx[i]<n){
					list.at(tmp_idx[i]) = n;
					tmp_idx[i] = n;
				}
			}
        }
	}
	list.at(ptr_out) = n;
	return 0;
}



void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i, unsigned int idx_out,
             const unsigned int ni, double *density_e, std::vector<double>& v_avg_x, std::vector<double>& v_avg_y,
             std::vector<double>& v_avg_z, std::vector<double>& v_rms_t, std::vector<double>& v_rms_l) {


    v_avg_x.clear();
    v_avg_x.resize(ni, 0);
    v_avg_y.clear();
    v_avg_y.resize(ni, 0);
    v_avg_z.clear();
    v_avg_z.resize(ni, 0);
    v_rms_l.clear();
    v_rms_l.resize(ni, 0);
    v_rms_t.clear();
    v_rms_t.resize(ni, 0);


    //Ions outside the electron bunch.
    while(list_i.at(idx_out)!=ni) {
        density_e[idx_out] = 0;
        v_avg_x.at(idx_out) = 0;
        v_avg_y.at(idx_out) = 0;
        v_avg_z.at(idx_out) = 0;
        v_rms_t.at(idx_out) = 0;
        v_rms_l.at(idx_out) = 0;
        idx_out = list_i.at(idx_out);
    }

    //Loop through the tree for ions inside the electron bunch.
    for (auto& box: tree) {
        if(box.n_ion>0) {
            double box_size = box.box_size;
            int n_e = box.n_ptcl;
            double d = n_e/(box_size*box_size*box_size); //density
            unsigned long int idx = box.first_ptcl;
            double vx_avg, vy_avg, vz_avg;
            vx_avg = 0;
            vy_avg = 0;
            vz_avg = 0;
            while(list_e.at(idx)!=ne) {
                vx_avg += vx.at(idx);
                vy_avg += vy.at(idx);
                vz_avg += vz.at(idx);
                idx = list_e.at(idx);
            }
            vx_avg /= n_e;
            vy_avg /= n_e;
            vz_avg /= n_e;

            idx = box.first_ptcl;
            double vx_rms, vy_rms, vz_rms;
            vx_rms = 0;
            vy_rms = 0;
            vz_rms = 0;
            while(list_e.at(idx)!=ne) {
                double dvx = vx.at(idx) - vx_avg;
                double dvy = vy.at(idx) - vy_avg;
                double dvz = vz.at(idx) - vz_avg;
                vx_rms += dvx*dvx;
                vy_rms += dvy*dvy;
                vz_rms += dvz*dvz;
                idx = list_e.at(idx);
            }

            vx_rms /= n_e;
            vy_rms /= n_e;
            vz_rms /= n_e;
            double vtr_rms = sqrt(vx_rms+vy_rms);
            vz_rms = sqrt(vz_rms);

            idx = box.first_ion;
            while(idx!=ni) {
                density_e[idx] = d;
                v_avg_x.at(idx) = vx_avg;
                v_avg_y.at(idx) = vy_avg;
                v_avg_z.at(idx) = vz_avg;
                v_rms_t.at(idx) = vtr_rms;
                v_rms_l.at(idx) = vz_rms;
                idx = list_i.at(idx);
            }
        }
    }
}


void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i, unsigned int idx_out,
             const unsigned int ni, double *density_e,std::vector<double>& v_avg_z, std::vector<double>& v_rms_t,
             std::vector<double>& v_rms_l) {

    v_avg_z.clear();
    v_avg_z.resize(ni, 0);
    v_rms_l.clear();
    v_rms_l.resize(ni, 0);
    v_rms_t.clear();
    v_rms_t.resize(ni, 0);


    //Ions outside the electron bunch.
    while(list_i.at(idx_out)!=ni) {
        density_e[idx_out] = 0;
        v_avg_z.at(idx_out) = 0;
        v_rms_t.at(idx_out) = 0;
        v_rms_l.at(idx_out) = 0;
        idx_out = list_i.at(idx_out);
    }

    //Loop through the tree for ions inside the electron bunch.
    for (auto& box: tree) {
        if(box.n_ion>0) {
            double box_size = box.box_size;
            int n_e = box.n_ptcl;
            double d = n_e/(box_size*box_size*box_size); //density
            unsigned long int idx = box.first_ptcl;
//            double vx_avg, vy_avg, vz_avg;
            double vz_avg;
            vz_avg = 0;
            while(list_e.at(idx)!=ne) {
                vz_avg += vz.at(idx);
                idx = list_e.at(idx);
            }
            vz_avg /= n_e;

            idx = box.first_ptcl;
            double vx_rms = 0, vy_rms = 0, vz_rms = 0;
            while(list_e.at(idx)!=ne) {
                double dvz = vz.at(idx) - vz_avg;
                vx_rms += vx.at(idx)*vx.at(idx);
                vy_rms += vy.at(idx)*vy.at(idx);
                vz_rms += dvz*dvz;
                idx = list_e.at(idx);
            }

            vx_rms /= n_e;
            vy_rms /= n_e;
            vz_rms /= n_e;
            double vtr_rms = sqrt(vx_rms+vy_rms);
            vz_rms = sqrt(vz_rms);

            idx = box.first_ion;
            while(idx!=ni) {
                density_e[idx] = d;
                v_avg_z.at(idx) = vz_avg;
                v_rms_t.at(idx) = vtr_rms;
                v_rms_l.at(idx) = vz_rms;
                idx = list_i.at(idx);
            }
        }
    }
}



void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i,
             unsigned int idx_out, const unsigned int ni, double *density_e,
             std::vector<double>& v_rms_t, std::vector<double>& v_rms_l) {

    v_rms_l.clear();
    v_rms_l.resize(ni, 0);
    v_rms_t.clear();
    v_rms_t.resize(ni, 0);

    //Ions outside the electron bunch.
    while(list_i[idx_out]!=ni) {
        density_e[idx_out] = 0;
        v_rms_t[idx_out] = 0;
        v_rms_l[idx_out] = 0;
        idx_out = list_i[idx_out];
    }

    //Loop through the tree for ions inside the electron bunch.
    for (auto& box: tree) {
        if(box.n_ion>0) {
            double box_size = box.box_size;
            int n_e = box.n_ptcl;
            double d = n_e/(box_size*box_size*box_size); //density
            unsigned long int idx = box.first_ptcl;
            double vx_rms, vy_rms, vz_rms;
            vx_rms = 0;
            vy_rms = 0;
            vz_rms = 0;
//            std::cout<<"======================="<<std::endl;
//            std::cout<<n_e<<std::endl;

            while(list_e.at(idx)!=ne) {
//                std::cout<<idx<<' '<<vx.at(idx)<<' '<<vy.at(idx)<<' '<<vz.at(idx)<<std::endl;
                vx_rms += vx.at(idx)*vx.at(idx);
                vy_rms += vy.at(idx)*vy.at(idx);
                vz_rms += vz.at(idx)*vz.at(idx);
                idx = list_e.at(idx);
            }

            vx_rms /= n_e;
            vy_rms /= n_e;
            vz_rms /= n_e;
            double vtr_rms = sqrt(vx_rms+vy_rms);
            vz_rms = sqrt(vz_rms);

            idx = box.first_ion;
            while(idx!=ni) {
                density_e[idx] = d;
                v_rms_t.at(idx) = vtr_rms;
                v_rms_l.at(idx) = vz_rms;
                idx = list_i.at(idx);
            }
        }
    }
}

//output a box for debug
std::ostream& operator<<(std::ostream& os, Box& box){
	os<<"Box center: ";
	for(int i=0; i<3; ++i) os<<box.center[i]<<' ';
	os<<endl;
	os<<"Box size: " << box.box_size <<endl;
	os<<"Parent box: "<< box.parent <<endl;
	os<<"Number of child boxes: "<< box.n_child << endl;
	os<<"Child boxes: ";
	for(int i=0; i<8; ++i) {
		if(box.child[i]>0) os<< box.child[i] <<' ';
	}
//	for(int i=0; i<box.n_child; ++i) {
//		os<< box.child[i] <<' ';
//	}
	os << endl;
	os<<"Particles: " << box.n_ptcl <<' ' << box.first_ptcl<<endl;
	os<<"Ions: " << box.n_ion <<' ' << box.first_ion<<endl;
	return os;
}

//output colleague for debug
std::ostream& operator<<(std::ostream& os, Colleague& clg){
	os<<clg.clg[0]<<' ';
	for(int i=1; clg.clg[i]>0; ++i) os<<clg.clg[i]<<' ';
	os<<endl;
	return os;
}

