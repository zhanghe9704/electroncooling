#ifndef ARBITRARY_ELECTRON_BEAM_H
#define ARBITRARY_ELECTRON_BEAM_H

#ifndef BOX_H
#define BOX_H

#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

typedef struct Box{
	double center[3] = {0,0,0};
	unsigned long int parent = 0;
	unsigned long int child[8] = {0};
	unsigned long int first_ptcl = 0;
	unsigned int n_ptcl = 0;
	int n_child = 0;
	double box_size = 0;
	unsigned long int first_ion = 0;
	unsigned int n_ion = 0;
} Box;

typedef struct Colleague{
	unsigned long int clg[28] = {0}; //At most 27 colleagues, the last 0 means the end.
} Colleague;

//unsigned long int load_electrons(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
//                                 std::vector<double>& vy, std::vector<double>& vz, unsigned long int n, std::string filename,
//                                int line_skip = 0);
unsigned long int load_electrons(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
                                 std::vector<double>& vy, std::vector<double>& vz, std::string filename, long int n,
                                int skip = 0, bool binary = false, int n_buffer = 1000);

//int create_e_tree(double * x, double * y, double * z, const unsigned long int n, const unsigned int s, vector<Box> &tree,
//                unsigned long int * list);
int create_e_tree(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const unsigned long int n,
                  const unsigned int s, vector<Box> &tree, std::vector<unsigned long int>& list);

//int create_ion_tree(double * x, double * y, double * z, const unsigned long int n, vector<Box> &tree,
//                unsigned long int * list, unsigned long int &idx_out);
int create_ion_tree(double * x, double * y, double * z, const unsigned int n, vector<Box> &tree,
                std::vector<unsigned int>& list, unsigned int &idx_out);

//void density(vector<Box> &tree, unsigned long int *list_e, double *vx, double *vy,
//             double *vz, const unsigned long int ne,  unsigned long int *list_i,
//             unsigned long int idx_out, const unsigned long int ni, double *density_e,
//             double *v_avg_x, double *v_avg_y, double *v_avg_z, double *v_rms_t, double *v_rms_l);
//
//void density(vector<Box> &tree, unsigned long int *list_e, double *vx, double *vy,
//             double *vz, const unsigned long int ne,  unsigned long int *list_i,
//             unsigned long int idx_out, const unsigned long int ni, double *density_e,
//             double *v_rms_t, double *v_rms_l);

void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i, unsigned int idx_out,
             const unsigned int ni, double *density_e, std::vector<double>& v_avg_x, std::vector<double>& v_avg_y,
             std::vector<double>& v_avg_z, std::vector<double>& v_rms_t, std::vector<double>& v_rms_l);

void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i, unsigned int idx_out,
             const unsigned int ni, double *density_e,std::vector<double>& v_avg_z, std::vector<double>& v_rms_t,
             std::vector<double>& v_rms_l);

void density(vector<Box> &tree, std::vector<unsigned long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const unsigned long int ne,  std::vector<unsigned int>& list_i,
             unsigned int idx_out, const unsigned int ni, double *density_e,
             std::vector<double>& v_rms_t, std::vector<double>& v_rms_l);

std::ostream& operator<<(std::ostream& os, Box& box);
std::ostream& operator<<(std::ostream& os, Colleague& clg);
#endif

#endif // ARBITRARY_ELECTRON_BEAM_H
