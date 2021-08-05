//---------------------------------------------------------------//
// GCAF= Guiding Center Analytic Field code by R. Saavedra       //
// Email of the author: saavestro@gmail.com                      //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git  //
//---------------------------------------------------------------//
#ifndef DECLAR_H
#define DECLAR_H

//CONSTANTS
const double pi= 3.14159265;      //pi number
const double keverg= 1.6022e-9;   //keV to erg
const double btkg= 1.0e+4;        //10 kG
const double rq= 4.8032e-10;      //elec charge
const double rc= 2.9979e+10;      //speed of light
const double mp= 1.6726e-24;      //proton mass
double sdum= 0.0;

//INPUT DATA VARIABLES, DEFAULTS ARE GIVEN
std::string indir= "./data.in";
std::string outdir= "./";
int    s_simu= 1;      //input simulation swtich traj/poin/ense
bool   s_conf= false;  //input configuration switch tkmk/stell
bool   s_coll= false;  //input collisions switch on/off
bool   s_pert= false;  //input perturbation switch on/off
int    i_part= 1;      //input number of particles
int    i_tran= 10;     //input number of toroidal transits
int    i_sdel= 1000;   //input step delta
bool   s_sfix= false;  //input step fix switch on/off
double d_dler= 1.0e-6; //input step fix en delta
int	   i_traj= 10000;  //input number of traj points
int    i_mome= 400;    //input number of disp points
int    i_dist= 4;      //input number of dist points
double r1= 0.4;        //poin lower bound
double r2= 0.6;        //poin upper bound
int    bm= 1;          //input on axis b field
int    rmaj= 150;      //input major radius
int    rmin= 50;       //input minor radius
double d_ampx= 1.0;    //input perturbation amplitude
int    nfrq= 0;        //input perturbation frequency
int    mmod= 2;        //input m mode number
int    nmod= 1;        //input n mode number
int    nphs= 0;        //input perturbation phase
double npto= 0.0;      //input pt in keV
int    s_potp= 0;      //input pt profile
int    mi= 1;          //ions mass in mp
int	   zi= 1;          //ions charge in protons charge
int    ntmpi= 1;       //ions temp in keV
int	   ntmpe= 1;       //elec temp in keV
bool   s_temp= 0;      //temp profile flat/gaus
double deno= 2.0e+14;  //plasma dens in cm^-3
bool   s_dens= 0;      //dens profile flat/para
int	   mt= 1;          //test mass in mp
int	   zt= 1;          //test charge in protons charge
int    ptchgrad= 0;    //input pitch in grades
double i_engo = 1.0;   //input initial en in keV
bool   s_ends= 0;      //input en distribution mon/max
int    nrado= rmin/2;  //input initial rad in cm
int    thtgrad= 0;     //theta in grades
int    zetgrad= 0;     //zeta in grades
int    cs= 1;          //test type ions/elec
std::string lngargtraj= "traj.plt    tim rad tht zet par per eng pch rcy zcy";
std::string lngargpoin= "poin.plt    rad tht rcy zcy";
std::string lngargmome= "disp.plt    tim pol eng pzt";
std::string lngargdist= "dist.plt    rad eng pch";
std::string lngargcoef= "coef.plt    den dfr dfe dfp sdr sde sdp";
std::array<int,5> lngargsize;

//SIMULATION
double trun;    //total simulation time
double omego;   //characteristic time freq
double leno;    //characteristic length
double enon;    //normalized initial en
double enaxo;   //on axis en
double dto;     //initial time step
std::array<double,MAX_PART> dt;  //time step
std::array<double,MAX_PART> tim; //time variable

//FIELD
double prtfrq;                //perturbation frequency
double pto= 0.0;              //pt magnitude
std::array<double,MAX_PART> bf;     //b field
std::array<double,MAX_PART> dbdp;   //b pol derivative
std::array<double,MAX_PART> dbdt;   //b tht derivative
std::array<double,MAX_PART> dbdz;   //b zet derivative
std::array<double,MAX_PART> al;     //alfa perturbation
std::array<double,MAX_PART> dadp;   //al pol derivative
std::array<double,MAX_PART> dadt;   //al tht derivative
std::array<double,MAX_PART> dadz;   //al zet derivative
std::array<double,MAX_PART> dadtim; //al zet derivative
std::array<double,MAX_PART> pt;     //electrostatic potential
std::array<double,MAX_PART> dptdp;  //pt pol derivative
std::array<double,MAX_PART> dptdt;  //pt tht derivative
std::array<double,MAX_PART> dptdz;  //pt zet derivative
std::array<double,MAX_PART> cg;     //poloidal current
std::array<double,MAX_PART> cgp;    //cg pol derivative
std::array<double,MAX_PART> ci;     //toroidal current
std::array<double,MAX_PART> cip;    //ci pol derivative
std::array<double,MAX_PART> qf;     //q function

//TEST PARTICLE
double ptcho;                 //initial pitch in cos
double eno;                   //initial en in keV
double polo;                  //initial poloidal flux
double theto;                 //initial theta
double zeto;                  //initial zeta
double polw;                  //on wall poloidal flux
std::array<double,MAX_PART> ptch;   //pitch
std::array<double,MAX_PART> en;     //energy
std::array<double,MAX_PART> pol;    //poloidal flux
std::array<double,MAX_PART> tht;    //theta
std::array<double,MAX_PART> zet;    //zeta
std::array<double,MAX_PART> rho;    //parallel gyroradius
std::array<double,MAX_PART> mu;     //b moment
std::array<double,MAX_PART> pz;     //zeta conjugate momenta
std::array<double,MAX_PART> toten;  //total en
std::array<double,MAX_PART> totenn; //total en normalized

//ENSEMBLE
//static const int dim_ens= 400;
std::array<double,MAX_MOMEPNTS> tim_moments; //time for disp
std::array<double,MAX_MOMEPNTS> rd1_av;  //rad av
std::array<double,MAX_MOMEPNTS> rd2_av;  //rad sq av
std::array<double,MAX_MOMEPNTS> rd3_av;  //rad cu av
std::array<double,MAX_MOMEPNTS> rd4_av;  //rad fo av
std::array<double,MAX_MOMEPNTS> en1_av;  //en av
std::array<double,MAX_MOMEPNTS> en2_av;  //en sq av
std::array<double,MAX_MOMEPNTS> en3_av;  //en cu av
std::array<double,MAX_MOMEPNTS> en4_av;  //en fo av
std::array<double,MAX_MOMEPNTS> pz1_av;  //pz av
std::array<double,MAX_MOMEPNTS> pz2_av;  //pz sq av
std::array<double,MAX_MOMEPNTS> pz3_av;  //pz cu av
std::array<double,MAX_MOMEPNTS> pz4_av;  //pz fo av
std::array<double,MAX_MOMEPNTS> rd1;     //rad one part
std::array<double,MAX_MOMEPNTS> rd2;     //rad sq one part
std::array<double,MAX_MOMEPNTS> rd3;     //rad cu one part
std::array<double,MAX_MOMEPNTS> rd4;     //rad fo one part
std::array<double,MAX_MOMEPNTS> en1;     //en one part
std::array<double,MAX_MOMEPNTS> en2;     //en sq one part
std::array<double,MAX_MOMEPNTS> en3;     //en cu one part
std::array<double,MAX_MOMEPNTS> en4;     //en fo one part
std::array<double,MAX_MOMEPNTS> pz1;     //pz one part
std::array<double,MAX_MOMEPNTS> pz2;     //pz sq one part
std::array<double,MAX_MOMEPNTS> pz3;     //pz cu one part
std::array<double,MAX_MOMEPNTS> pz4;     //pz fo one part
std::array<double,MAX_MOMEPNTS> rprv;    //rad previous
std::array<double,MAX_PART> dist_tim;
std::array<std::array<double,MAX_DISTPNTS>,MAX_PART> dist_rad;
std::array<std::array<double,MAX_DISTPNTS>,MAX_PART> dist_eng;
std::array<std::array<double,MAX_DISTPNTS>,MAX_PART> dist_pch;

//RECORDS
int   flux;                                //particle flux
double zplo;                               //poin zeta plane
std::array<double,MAX_PART> pold;          //old pol for poin
std::array<double,MAX_PART> told;          //old tht for poin
std::array<double,MAX_PART> zold;          //old zet for poin
std::array<double,RCRDTRAJ_SIZE> rcrdtraj; //record traj points
std::vector<double> rcrdpoin;              //record traj points
double* pnt_traj[MAX_DATACOLS];
double* tgt_traj[MAX_DATACOLS];
double* pnt_poin[MAX_DATACOLS];
double* tgt_poin[MAX_DATACOLS];
double adress[MAX_DATACOLS];
std::fstream traj;
std::fstream poin;
std::fstream mome;
std::fstream dist;
std::array<std::fstream,MAX_DISTPNTS> dist_n;
std::fstream dati;
std::fstream dran;
std::fstream dlog;

#endif
