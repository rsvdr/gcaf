//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#ifndef DEFINE_H
#define DEFINE_H

// Input data
std::string outdir= "";                        // output directory
struct Switches swt;                           // switches variables
struct Simulation sim;                         // simulation variables
struct Field field;                            // field variables
struct Species species[MAX_SPECIES];           // impurities species variables
struct Species plasma;                         // plasma variables
int total_species= 1;                          // total number of species
struct Particle test;                          // test particle variables
struct Code code;                              // code variables

// Field array variables
std::array<double,MAX_PART> bf;                // b field
std::array<double,MAX_PART> bf_pol;            // b pol derivative
std::array<double,MAX_PART> bf_tht;            // b tht derivative
std::array<double,MAX_PART> bf_zet;            // b zet derivative
std::array<double,MAX_PART> al;                // alfa perturbation
std::array<double,MAX_PART> al_pol;            // al pol derivative
std::array<double,MAX_PART> al_tht;            // al tht derivative
std::array<double,MAX_PART> al_zet;            // al zet derivative
std::array<double,MAX_PART> al_tim;            // al zet derivative
std::array<double,MAX_PART> pt;                // electrostatic potential
std::array<double,MAX_PART> pt_pol;            // pt pol derivative
std::array<double,MAX_PART> pt_tht;            // pt tht derivative
std::array<double,MAX_PART> pt_zet;            // pt zet derivative
std::array<double,MAX_PART> cg;                // poloidal current
std::array<double,MAX_PART> cg_pol;            // cg pol derivative
std::array<double,MAX_PART> ci;                // toroidal current
std::array<double,MAX_PART> ci_pol;            // ci pol derivative
std::array<double,MAX_PART> qf;                // q profile
std::array<double,MAX_PART> qf_pol;            // q profile derivative

// Integration array variables
std::array<double,MAX_PART> dt;                // time step
std::array<double,MAX_PART> tim;               // time variab
std::array<double,MAX_PART> ptch;              // pitch
std::array<double,MAX_PART> en;                // energy
std::array<double,MAX_PART> pol;               // poloidal flux
std::array<double,MAX_PART> tht;               // theta
std::array<double,MAX_PART> zet;               // zeta
std::array<double,MAX_PART> rho;               // parallel gyroradius
std::array<double,MAX_PART> mu;                // b moment
std::array<double,MAX_PART> pz;                // zeta conjugate momenta
std::array<double,MAX_PART> toten;             // total

// Ensemble array variables
int total_slices= 0;                           // number of slices
int ip_size= 0;                                // interpolated data size
std::array<double,MAX_POINPNTS> rad_ip;        // rad interpolated data
std::array<double,MAX_POINPNTS> tht_ip;        // tht interpolated data
std::array<int,MAX_PART>        flux;          // particle flux
std::array<double,MAX_PART>     rad_prev;      // rad previous

#endif
