//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#ifndef SHARED_H
#define SHARED_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <exception>


// Bounds
static const std::size_t MAX_DATACOLS= 1e1; //max number of data columns
static const std::size_t MAX_TRAJPNTS= 1e4; //max traj points
static const std::size_t MAX_MOMEPNTS= 1e4; //max mome points
static const std::size_t MAX_POINPNTS= 1e6; //max mome points
static const std::size_t MAX_DISTPNTS= 4e1; //max dist points
static const std::size_t MAX_PART    = 1e4; //max number of particles
static const std::size_t MAX_SDEL    = 4e3; //max step delta
static const std::size_t MAX_TRAN    = 2e3; //max number of transits
static const std::size_t MAX_SPECIES = 1e1; //max number of species
static const std::size_t RCRDTRAJ_SIZE=     //max rcrd_traj variable size
                            MAX_DATACOLS*MAX_TRAJPNTS;
static const std::size_t RCRDPOIN_SIZE=     //max rcrd_poin variable size
                            MAX_DATACOLS*MAX_POINPNTS;
static const std::size_t RCRDMOME_SIZE=      //max rcrd_mome variable size
                            MAX_DATACOLS*MAX_MOMEPNTS;

// Constants
static const double KEV_TO_ERG   = 1.60217656e-09; // keV to erg
static const double TSLA_TO_GAUS = 1.00000000e+04; // 10 kG
static const double CSPEED       = 2.99792458e+10; // speed of light
static const double CHRG_PROT    = 4.80320451e-10; // proton charge
static const double MASS_PROT    = 1.67262192e-24; // proton mass
static const double MASS_ELEC    = 9.10938370e-28; // electron mass

// Types
enum sim_type_t // simulation types
    {SIM_TRAJ,SIM_ENSE,SIM_POIN};
enum dev_type_t // device types
    {DEV_TKMK,DEV_STEL};
enum pot_prof_t // potential profile types
    {POT_CONS,POT_PRES};
enum tem_prof_t // temperature profile types
    {TEM_FLAT,TEM_GAUS};
enum den_prof_t // density profile types
    {DEN_FLAT,DEN_PARA};
enum chg_type_t // charge profile types
    {CHG_IONS= 1,CHG_ELEC= -1};
enum eng_dist_t // energy distribution types
    {ENG_MONO,ENG_MAXW};
enum pch_dist_t // pitch angle distribution types
    {PCH_PAR,PCH_ANT,PCH_ISO,PCH_FIX};

extern std::string outdir;

struct Switches{ // Input switches
    sim_type_t simu= SIM_TRAJ;       //sim type traj/poin/ense
    dev_type_t conf= DEV_TKMK;       //conf tkmk/stell
    bool       coll= false;          //collisions
    bool       pert= false;          //perturbation
    bool       scrn= false;          //screen log
    bool       verb= false;          //verbose
    bool       flux= false;          //fluxxz
};
extern struct Switches swt;

struct Simulation{ // Input simulation variables
    int        i_part= 1;            //number of particles
    int        i_tran= 10;           //number of toroidal transits
    int        i_sdel= 1000;         //step delta
    int	       i_traj= 1000;         //number of traj points
    int        i_mome= 400;          //number of disp points
    int        i_dist= 4;            //number of dist points
    double     d_rlow= 0.4;          //poin lower bound
    double     d_rupp= 0.6;          //poin upper bound
    int        i_surf= 0;            //const zeta plane in grades
};
extern struct Simulation sim;

struct Field{ // Input field variables
    int        i_baxi= 1;            // on axis field magnitude
    int        i_qaxi= 1;            // on axis q profile
    int        i_qwll= 5;            // on wall q profile
    int        i_rmaj= 150;          // major radius
    int        i_rmin= 50;           // minor radius
    double     d_ampx= 1.0e-4;       // perturbation amplitude
    int        i_freq= 0;            // perturbation frequency
    int        i_phas= 0;            // perturbation phase
    int        i_mmod= 2;            // m mode number
    int        i_nmod= 1;            // n mode number
    double     d_epot= 0.0;          // pt in keV
    pot_prof_t potp= POT_CONS;       // pt profile
};
extern struct Field field;

struct Species{ // Input impurities species variables
    chg_type_t sign= CHG_IONS;       // charge type
    int        i_mass= 1;            // mass in proton mass
    int        i_chrg= 1;            // charge in proton charge
    tem_prof_t temp_prof= TEM_FLAT;  // temperature distribution
    int        i_temp= 1;            // temperature in keV
    den_prof_t dens_prof= DEN_FLAT;  // density distribution
    double     d_dens= 2.0e+14;      // density in cm^-3
};
extern struct Species species[MAX_SPECIES];
extern struct Species plasma;
extern int total_species;


struct Particle{ // Input test particle variables
    chg_type_t sign= CHG_IONS;       // charge type
    int        i_mass= 1;            // mass in proton mass
    int        i_chrg= 1;            // charge in proton charge
    double     d_radi= 0.5;          // radial position in rmin
    int        i_thet= 0;            // theta in grades
    int        i_zeta= 0;            // zeta in grades
    double     d_ener= 1.0;          // energy in keV
    eng_dist_t ener_dist= ENG_MONO;  // energy distribution
    int        i_ptch= 0;            // pitch angle in grades
    pch_dist_t ptch_dist= PCH_PAR;   // pitch angle distribution
};
extern struct Particle test;

struct Code{ // Normalized variables
    double     d_trun= 0.0;          // total simulation time
    double     d_freq= 0.0;          // characteristic time freq
    double     d_leng= 0.0;          // characteristic length
    double     d_enax= 0.0;          // characteristic energy
    double     d_dt= 0.0;            // initial time step
    double     d_pert_freq= 0.0;     // perturbation frequency
    double     d_ptch= 0.0;          // pt magnitude
    double     d_en= 0.0;            // initial pitch in cos
    double     d_poli= 0.0;          // initial en in keV
    double     d_thet= 0.0;          // initial poloidal flux
    double     d_zeta= 0.0;          // initial theta
    double     d_polw= 0.0;          // initial zeta
    double     d_pt= 0.0;            // on wall poloidal flux
    double     d_zeta_poin= 0.0;     // poin zeta plane
    double     d_test_mass= 0.0;     // test particle mass
    double     d_test_chrg= 0.0;     // test particle charge
};
extern struct Code code;

// Field variables
extern std::array<double,MAX_PART> bf;
extern std::array<double,MAX_PART> bf_pol;
extern std::array<double,MAX_PART> bf_tht;
extern std::array<double,MAX_PART> bf_zet;
extern std::array<double,MAX_PART> al;
extern std::array<double,MAX_PART> al_pol;
extern std::array<double,MAX_PART> al_tht;
extern std::array<double,MAX_PART> al_zet;
extern std::array<double,MAX_PART> al_tim;
extern std::array<double,MAX_PART> pt;
extern std::array<double,MAX_PART> pt_pol;
extern std::array<double,MAX_PART> pt_tht;
extern std::array<double,MAX_PART> pt_zet;
extern std::array<double,MAX_PART> cg;
extern std::array<double,MAX_PART> cg_pol;
extern std::array<double,MAX_PART> ci;
extern std::array<double,MAX_PART> ci_pol;
extern std::array<double,MAX_PART> qf;
extern std::array<double,MAX_PART> qf_pol;

// Integration variables
extern std::array<double,MAX_PART> dt;
extern std::array<double,MAX_PART> tim;
extern std::array<double,MAX_PART> ptch;
extern std::array<double,MAX_PART> en;
extern std::array<double,MAX_PART> pol;
extern std::array<double,MAX_PART> tht;
extern std::array<double,MAX_PART> zet;
extern std::array<double,MAX_PART> rho;
extern std::array<double,MAX_PART> mu;
extern std::array<double,MAX_PART> pz;
extern std::array<double,MAX_PART> toten;

// Ensemble variables
extern int                             total_slices;
extern int                             ip_size;
extern std::array<double,MAX_POINPNTS> rad_ip;
extern std::array<double,MAX_POINPNTS> tht_ip;
extern std::array<int,MAX_PART>        flux;
extern std::array<double,MAX_PART>     rad_prev;

#endif
