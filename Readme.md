#Introduction of the new electron cooling simulation code#

## Environment  ##
This code is devloped by C\++. Users need to write their own "main" program. To compile the code, please make sure your compiler supports c\++ 11 and  the option for c\++ 11 is turned on. For example, when using GCC, one needs to add "-std=c\++11" into the compiler options. A Code Blocks Project (cbp) file is included in the subfolder cooling_release, and the Code Blocks IDE with a proper compiler can be used to compiler the project. Makefiles for multiple platform can be generated using the tool "cbp2make". 

## Constants ##
The following physical and mathematical constants are defined in "__constants.h__":
~~~~c++
const double k_ke = 8.9875517873681764E9;   //Coulomb's constant, in N*m^2/C^2
const double k_c = 299792458.0;             //speed of light, in m/s
const double k_e = 1.602176565E-19;         //Elementary charge, in C
const double k_pi = 3.1415926535897932384626;
const double k_u = 931.49406121;            //Atomic mass unit, in MeV/c^2
const double k_me = 0.510998928;            //electron mass, in MeV/c^2
~~~~

## Global variables ##

 For intrabeam scattering (IBS) and/or electron cooling dynamic simulations, it is required to use the following variables to save configurations of respective calculations. (Will be explained in the following sections.) The declaration of them, as shown in the code block below, should be included in the main program.  For rate calculation, they can be used, but not required. 

~~~~c++
extern IBSParas * ibs_paras;
extern EcoolRateParas * ecool_paras;
extern ForceParas * force_paras;
extern DynamicParas * dynamic_paras;
//Only used for IBS simulation (without electron cooling) using model beam method
extern Twiss * twiss_ref;	
extern int n_ion_model;		
~~~~

## Define the ion beam ##
The class __Beam__ is provided in the "__beam.h__". An ion beam can be defined using this class. The user should set value for the following paramters: 
* charge number of the ion 
* mass number of the ion 
* kinetic energy in MeV
* normalized horizontal emittance in m*rad
* normalized vertical emittance in m*rad
* mommentum spread dp/p
* rms bunch length (for bunched beam only)
* number of ions in the beam (for coasting beam) or in each bunch (for bunched beam)

The following code defines a bunched proton beam:
~~~~c++
double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl, A;
int Z;
Z = 1;				//Charge number is 1. 
m0 = 938.272;		//Mass in MeV
A = m0/k_u			//Mass number
KE = 250e3;			//Kinetic energy in MeV
emit_nx0 = 1e-6;	//Transverse normalized emittance in m*rand
emit_ny0 = 0.5e-6;	//Vertical normalized emittance in m*rad
dp_p0 = 0.0007;		//Momentum spread dp/p
N_ptcl = 6.56E9;	//Number of protons per bunch
sigma_s0 = 2E-2;	//RMS bunch length
//Define the proton beam as p_beam
Beam p_beam(Z, A, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);
~~~~
A coasting proton beam can be defined in the same way without the rms bunch length included, such as 
~~~~c++
Beam p_beam(Z, A, KE, emit_nx0, emit_ny0, dp_p0, N_ptcl);
~~~~



## Define the ring ##
The class __Lattice__ and the class __Ring__ are provided in the "__ring.h__". 

The TWISS parameters at different position of the machine in MADX tsf format can be read into the class Lattice. The TWISS parameters should be saved in the following sequency: "s, bet_x, alf_x, mu_x, d_x, dp_x, bet_y, alf_y, mu_y, d_y, dp_y". The following code define a class lattice that saves the TWISS parameters:
~~~~c++
//The TWISS parameters are saved in MEICColliderRedesign1IP.tfs
std::string filename = "MEICColliderRedesign1IP.tfs";  
Lattice lattice(filename);	//Load the TWISS parameters to lattice
~~~~
The ring can be defined once the TWISS parameters and the ion beam are ready. 
~~~~c++
Ring ring(lattice, p_beam);
~~~~



## Intrabeam scattering (IBS) rate calculation

The IBS rate can be calculated, after the ring and the ion beam are define. The IBS rate is calculated using Martini model. One needs to set up the grid number for the 3D integral in Martini formula. 

~~~~c++
int nu = 100;
int nv = 100;
int nz = 40;
ibs_paras = new IBSParas(nu, nv, nz);
~~~~

Coulomb logarithm can be used to replace the longitudinal integral. 

~~~~c++
int nu = 100;
int nv = 100;
double log_c = 19.4;
ibs_paras = new IBSParas(nu, nv, log_c);
~~~~

The transverse coupling rate can be set between 0 to 1. 

~~~~c++
ibs_paras->set_k(0.2); //Set the transverse coupling rate to be 0.2. Default is 0. 
~~~~

Now the IBS rate can be calculated as follows. 

~~~~c++
double rx_ibs, ry_ibs, rz_ibs;	//Variables to save the IBS rate results
config_ibs(lattice);
ibs_rate(lattice, p_beam, *ibs_paras, rx_ibs, ry_ibs, rz_ibs);
end_ibs();
//Print the IBS rate to screen. 
std::cout<<"ibs rate: "<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
~~~~



## Define the cooler ##

The __Cooler__ class is provide in "cooler.h". The following parameters are required to define an electron cooler: the length of the cooler, the number of coolers, the magnetic field inside the cooler, and the transverse beta functions at the cooler. One can also specifies the dispersion functions, the alpha functions, and the derivative of the dispersion functions at the cooler.  If not set, the default values for them are zeros. 

~~~~c++
double cooler_length = 10;
double n_section = 1;
double magnetic_field = 0.1;
double beta_h = 10;
double beta_v = 10;
double dis_h = 0;	//Dispersion functions at the cooler can be set but not required.
double dis_v = 0;	
Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);
~~~~



## Define the electron beam ##

The __EBeamShape__ class in "beam.h" provides an interface for different shapes of the electron beam. Currently the following three shapes have been defined: coasting electron beam with uniform density and cylinder shape, Gaussian bunch, and uniform bunch. (Other shapes can be defined easily as a derived class of __EBeamShape__ too.)  The respective classes for them  are __UniformCylinder__ , __GaussianBunch__ and __UniformBunch__ , all of which are derived classes of __EBeamShape__ .  Users can set values of different parameters for each shape, as follows.

~~~~c++
//Define a coasting electron beam with uniform density and cylinder shape
double current = 2;		//Current of the electron beam in A
double radius = 0.008;	//Transverse radius of the electron beam
UniformCylinder uniform_cylinder(current, radius);
~~~~

~~~~c++
//Define a bunched electron beam with Gaussian distribution
double ne = 1e8; 		//Number of electrons in the bunch
double sigma_x = 1.5e-2;//RMS bunch size in the three dimensions in meters
double sigma_y = 1.5e-2;
double sigma_s = 2e-2;
GaussianBunch gaussian_bunch(ne, sigma_x, sigma_y, sigma_s);
~~~~

~~~~c++
//Define a bunched electron beam with uniform density and beer can shape
double current = 2;		//Current of the electron beam in A
double radius = 0.008;	//Transverse radius of the electron beam
double length = 2e-2	//Full length of the bunch in meters
UniformBunch uniform_bunch(current, radius, length);
~~~~

Then the electron beam in the cooler can be defined using the  __EBeam__ class,  provided in "beam.h" Users can set value of the temperature and the Lorentz factor $\gamma$  for the electron beam.  

~~~~c++
double gamma_e = p_beam.gamma();	//Electrons have the same velocity with the ions
double tmp_tr = 0.1;				//Transverse temperature
double tmp_long = 0.1;				//Longitudianl temperature
EBeam e_beam(gamma_e, tmp_tr, tmp_long, uniform_cylinder); //Coasting electron beam
~~~~



## Choose the friction force formula ##

## Electron cooling rate calculation ##

## Simulation of the electron cooling process ##

## Advanced topics ##

## Sample code ##

~~~~c++

~~~~

