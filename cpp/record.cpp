//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://github.com/Saavestro/gcafpp //
//----------------------------------------------------------------//
#include "shared.h"

struct Moments{ // Distribution moments variables
    std::array<double,MAX_MOMEPNTS> pow1;    // x^1
    std::array<double,MAX_MOMEPNTS> pow2;    // x^2
    std::array<double,MAX_MOMEPNTS> pow3;    // x^3
    std::array<double,MAX_MOMEPNTS> pow4;    // x^4
    std::array<double,MAX_MOMEPNTS> pow1_av; // x^1 average
    std::array<double,MAX_MOMEPNTS> pow2_av; // x^2 average
    std::array<double,MAX_MOMEPNTS> pow3_av; // x^3 average
    std::array<double,MAX_MOMEPNTS> pow4_av; // x^4 average
};
struct Moments mome_rd; // radial position moments
struct Moments mome_en; // energy moments
struct Moments mome_pz; // conjugate momenta
struct Distribution{ // Distribution variables
    std::array<double,MAX_PART>                            tim; // time
    std::array<std::array<double,MAX_DISTPNTS+1>,MAX_PART> rad; // radial pos
    std::array<std::array<double,MAX_DISTPNTS+1>,MAX_PART> tht; // theta
    std::array<std::array<double,MAX_DISTPNTS+1>,MAX_PART> zet; // zeta
    std::array<std::array<double,MAX_DISTPNTS+1>,MAX_PART> eng; // energy
    std::array<std::array<double,MAX_DISTPNTS+1>,MAX_PART> pch; // pitch angle
};
struct Distribution dist; // dsitribution

// Ensemble variables
std::array<double,MAX_MOMEPNTS>  mome_tim;

// Record variables
std::array<double,MAX_PART>      pol_old;   // old pol for poin
std::array<double,MAX_PART>      tht_old;   // old tht for poin
std::array<double,MAX_PART>      zet_old;   // old zet for poin
std::array<double,RCRDTRAJ_SIZE> rcrd_traj; // record traj points
std::array<double,RCRDMOME_SIZE> rcrd_flux; // record mome points
std::array<double,RCRDMOME_SIZE> rcrd_mome; // record flux points
std::vector<double>              rcrd_poin; // record poin points

// Function declarations
const double xproj(const double&, const double&, const double&);
const double yproj(const double&, const double&, const double&);
const double zproj(const double&, const double&, const double&);


void record_traj(const int& k,const int& skip){ //record traj points
	const double tim_dum= // time in s
		tim[k]/code.d_freq;
	const double rad_dum= // radial position in cm
		sqrt(2*pol[k])*code.d_leng;
	double tht_dum= tht[k];
	{ //bound -pi< tht< pi
		tht_dum= fmod(tht_dum,2*M_PI);
		if(tht_dum<-M_PI) tht_dum+= +2*M_PI;
		if(tht_dum>M_PI) tht_dum+= -2*M_PI;
	}
	double zet_dum= zet[k];
	{ //bound 0< zet< 2*pi
		zet_dum= fmod(zet_dum,2*M_PI);
		if(zet_dum<0) zet_dum+= +2*M_PI;
	}

	const double pparl= // parallel moment
		test.i_mass*code.d_leng*code.d_freq*rho[k]*bf[k]/CSPEED;
	double pperp= // perpendicular moment
		test.i_mass*code.d_leng*code.d_freq
		*sqrt(en[k] -0.5*rho[k]*rho[k]*bf[k]*bf[k])/CSPEED;
	if(pperp!=pperp) pperp= 0.0;
	const double xdum= // cartesian x
		xproj(pol[k],tht[k],zet[k]);
	const double ydum= // cartesian y
	 	yproj(pol[k],tht[k],zet[k]);
	const double zdum= // cartesian z
		zproj(pol[k],tht[k],zet[k]);
	const double bigr= // cylindrical R
	 	sqrt(xdum*xdum +ydum*ydum);

	rcrd_traj[0*sim.i_traj +skip]= tim_dum;
	rcrd_traj[1*sim.i_traj +skip]= rad_dum/field.i_rmin;
	rcrd_traj[2*sim.i_traj +skip]= tht_dum;
	rcrd_traj[3*sim.i_traj +skip]= zet_dum;
	// rcrd_traj[4*sim.i_traj +skip]= pparl;
	// rcrd_traj[5*sim.i_traj +skip]= pperp;
	rcrd_traj[4*sim.i_traj +skip]= en[k]*code.d_enax/KEV_TO_ERG;
	rcrd_traj[5*sim.i_traj +skip]= ptch[k];
	rcrd_traj[6*sim.i_traj +skip]= xdum;
	rcrd_traj[7*sim.i_traj +skip]= ydum;
	rcrd_traj[8*sim.i_traj +skip]= zdum;
}

void write_traj(){ // write traj.plt
	std::fstream traj;
	std::ofstream(outdir+"/traj.plt");
	traj.open(outdir+"/traj.plt",std::ios::app);
	const int numdat= 9; //number of data columns
	traj <<"#traj.plt    tim rad tht zet eng pch x y z\n";
	traj <<std::setprecision(8) <<std::scientific;
	for(int i=1;i<sim.i_traj;i++){
		for(int j=0;j<numdat;j++){
			traj << rcrd_traj[j*sim.i_traj +i] << "  ";
		}
		traj << "\n";
	}
	traj.close();
}

void record_flux(const int& k, const int& skip){ // record flux
	const double rad_dum= // radial position in cm
		code.d_leng*sqrt(2*pol[k]);
	double tht_dum= tht[k];
	{ //bound -pi< tht< pi
		tht_dum= fmod(tht_dum,2*M_PI);
		if(tht_dum<-M_PI) tht_dum+= +2*M_PI;
		if(tht_dum>M_PI) tht_dum+= -2*M_PI;
	}
	{ //bound 0< tht< ip_size
		tht_dum+= +M_PI;
		tht_dum*= ip_size/2/M_PI;
		if(tht_dum<0) tht_dum+= +ip_size;
		if(tht_dum>ip_size) tht_dum+= -ip_size;
	}
	double zet_dum= zet[k];
	{ //bound 0< zet< 2*pi
		zet_dum= fmod(zet_dum,2*M_PI);
		if(zet_dum<0) zet_dum+= +2*M_PI;
	}
	{ //bound 0< zet< 8
		zet_dum*= total_slices/2/M_PI;
		if(zet_dum<0) zet_dum+= +total_slices;
		if(zet_dum>total_slices) zet_dum+= -total_slices;
	}
	const int int_tht_dum= (int)tht_dum;
	const int int_zet_dum= (int)zet_dum;
	const double plane= rad_ip[int_zet_dum +int_tht_dum];
	if(rad_prev[k]/field.i_rmin < plane && rad_dum/field.i_rmin > plane)
		flux[k]+= 1;
	if(rad_prev[k]/field.i_rmin > plane && rad_dum/field.i_rmin < plane)
		flux[k]-= 1;
}


void write_flux(){ // write flux
	std::fstream flxf;
	std::ofstream(outdir+"/flux.plt");
	flxf.open(outdir+"/flux.plt",std::ios::app);
	for(int i=0;i<sim.i_mome;i++){
		flxf <<std::scientific <<std::setprecision(8)
			<< mome_tim[i]/code.d_freq <<"  " ;
		flxf <<std::scientific <<std::setprecision(8)
			<<1.0*rcrd_flux[i]/sim.i_part <<'\n';
	}
	flxf.close();
}

void record_mome(const int& k,const int& skip){ // record moments
	const double rad_dum= // radial position in cm
		code.d_leng*sqrt(2*pol[k]);
	mome_tim[skip]= tim[k];
	mome_rd.pow1[skip]= rad_dum;
	mome_rd.pow2[skip]= rad_dum*rad_dum;
	mome_rd.pow3[skip]= rad_dum*rad_dum*rad_dum;
	mome_rd.pow4[skip]= rad_dum*rad_dum*rad_dum*rad_dum;
	mome_en.pow1[skip]= en[k];
	mome_en.pow2[skip]= en[k]*en[k];
	mome_en.pow3[skip]= en[k]*en[k]*en[k];
	mome_en.pow4[skip]= en[k]*en[k]*en[k]*en[k];
	mome_pz.pow1[skip]= pz[k];
	mome_pz.pow2[skip]= pz[k]*pz[k];
	mome_pz.pow3[skip]= pz[k]*pz[k]*pz[k];
	mome_pz.pow4[skip]= pz[k]*pz[k]*pz[k]*pz[k];

	mome_rd.pow1_av[skip]+= mome_rd.pow1[skip];
	mome_rd.pow2_av[skip]+= mome_rd.pow2[skip];
	mome_rd.pow3_av[skip]+= mome_rd.pow3[skip];
	mome_rd.pow4_av[skip]+= mome_rd.pow4[skip];
	mome_en.pow1_av[skip]+= mome_en.pow1[skip];
	mome_en.pow2_av[skip]+= mome_en.pow2[skip];
	mome_en.pow3_av[skip]+= mome_en.pow3[skip];
	mome_en.pow4_av[skip]+= mome_en.pow4[skip];
	mome_pz.pow1_av[skip]+= mome_pz.pow1[skip];
	mome_pz.pow2_av[skip]+= mome_pz.pow2[skip];
	mome_pz.pow3_av[skip]+= mome_pz.pow3[skip];
	mome_pz.pow4_av[skip]+= mome_pz.pow4[skip];

	rcrd_flux[skip]+= flux[k];
}


const std::array<double,3> moments(
	const std::array<double,4>& mom){ // write moments helper function
	double av= mom[0]/sim.i_part;
	double sqav= mom[1]/sim.i_part;
	double cuav= mom[2]/sim.i_part;
	double foav= mom[3]/sim.i_part;
	double sdev= sqrt(fabs(sqav- av*av));
	return
    {
        sqav- av*av,  //variance
		(cuav -3*av*sqav +2*av*av*av)
        /sdev/sdev/sdev, //skewness
		(foav -4*av*cuav +6*sqav*av*av -3*av*av*av*av)
		/sdev/sdev/sdev/sdev //flatness
    };
}

void write_mome(){ // write moments
	std::fstream mome;
	std::ofstream(outdir+"/mome.plt");
	mome.open(outdir+"/mome.plt",std::ios::app);
	mome <<"#mome.plt    tim pol eng pzat\n";
	for(int i=sim.i_mome/10;i<sim.i_mome;i++){
		const std::array<double,3> rd_mome=
    		moments({mome_rd.pow1_av[i],
    				 mome_rd.pow2_av[i],
    				 mome_rd.pow3_av[i],
    				 mome_rd.pow4_av[i]});
		const std::array<double,3> en_mome=
    		moments({mome_en.pow1_av[i],
    			     mome_en.pow2_av[i],
    			     mome_en.pow3_av[i],
    			     mome_en.pow4_av[i]});
		const std::array<double,3> pz_mome=
    		moments({mome_pz.pow1_av[i],
    			     mome_pz.pow2_av[i],
    			     mome_pz.pow3_av[i],
    			     mome_pz.pow4_av[i]});

		 mome <<std::scientific <<std::setprecision(8)
		 	  <<mome_tim[i]/code.d_freq <<"  ";
		 for(int i=0;i<rd_mome.size();i++){
			 mome <<rd_mome[i] <<"  "
			      <<en_mome[i] <<"  "
				  <<pz_mome[i] <<"  ";
	     }
		 mome <<'\n';
	}
	mome.close();
}


void record_dist(const int& k, const int& skip){ //record distribution
	const double rad_dum= sqrt(2*pol[k])*code.d_leng/field.i_rmin;
	double tht_dum= tht[k];
	double zet_dum= zet[k];
	const double xdum= xproj(pol[k],tht[k],zet[k]);
	const double ydum= yproj(pol[k],tht[k],zet[k]);
	const double zdum= zproj(pol[k],tht[k],zet[k]);

	{//bound -pi< thet< pi
		tht_dum= fmod(tht_dum,2*M_PI);
		if(tht_dum<-M_PI) tht_dum+= +2*M_PI;
		if(tht_dum>M_PI) tht_dum+= -2*M_PI;
	}

	{ //bound 0< zet< 2*pi
		zet_dum= fmod(zet_dum,2*M_PI);
		if(zet_dum<0) zet_dum+= +2*M_PI;
	}

	dist.tim[k]= tim[k]/code.d_freq;
	dist.rad[k][skip]= rad_dum;
	dist.tht[k][skip]= tht_dum;
	dist.zet[k][skip]= zet_dum;
	// dist.rad[k][skip]= xdum;
	// dist.tht[k][skip]= ydum;
	// dist.zet[k][skip]= zdum;
	dist.eng[k][skip]= code.d_enax*en[k]/KEV_TO_ERG;
	dist.pch[k][skip]= ptch[k];
}


void write_dist(){ //write distribution
	std::array<std::fstream,MAX_DISTPNTS+1>dist_n;
	for(int i=0;i<sim.i_dist+1;i++){ // open sim.i_dist files + initial dist
		std::ofstream(outdir+"/dist_"+std::to_string(i)+".plt");
		dist_n[i].open(outdir+"/dist_"+std::to_string(i)+".plt",std::ios::app);
	}
	for(int i=0;i<sim.i_dist+1;i++){
		dist_n[i] << "dist_"+std::to_string(i)+".plt    tim" << '\n';
		dist_n[i] << dist.tim[i] <<'\n';
		dist_n[i] <<"# dist.plt    rad tht zet eng pch\n";
		for(int j=0;j<sim.i_part;j++)
			{dist_n[i] <<std::scientific <<std::setprecision(8)
					   <<dist.rad[j][i] <<"  "
					   <<dist.tht[j][i] <<"  "
					   <<dist.zet[j][i] <<"  "
				       <<dist.eng[j][i] <<"  "
				       <<dist.pch[j][i] <<'\n';}
	}
	for(int i=0;i<MAX_DISTPNTS+1;i++) dist_n[i].close();
}

void record_poin(const int& k){ // record poincare
	double pol_dum;
	double tht_dum;
	double zet_dum= (field.i_nmod*zet[k] -code.d_pert_freq*tim[k]);
	const double del= M_PI/1e3; //zeta plane infinitesimal width
	const double zpl= code.d_zeta_poin;     //zeta constant plane 0 < zpl <2*pi

	{ //bound 0< zet< 2*pi
		zet_dum= fmod(zet_dum,2*M_PI);
		if(zet_dum<0) zet_dum+= +2*M_PI;
	}

	{//check if particle has crossed zeta plane
		const double dum= zet_old[k];
	    if(dum>zpl+del || dum<zpl)
			{pol_old[k]= pol[k]; tht_old[k]= tht[k]; zet_old[k]= zet_dum; return;}
	}

	{ //linear interpolate and evaluate at plane
		pol_dum= (pol_old[k]*(zet_dum -zpl) +pol[k]*(zpl -zet_old[k]))
			/(zet_dum -zet_old[k]);
		tht_dum= (tht_old[k]*(zet_dum -zpl) +tht[k]*(zpl -zet_old[k]))
			/(zet_dum -zet_old[k]);
	}
	const double rad_dum= sqrt(2*pol[k])*code.d_leng/field.i_rmin;

	const double x_proj= xproj(rad_dum,tht_dum,zpl);
	const double y_proj= yproj(rad_dum,tht_dum,zpl);
	const double z_proj= zproj(rad_dum,tht_dum,zpl);
	const double r_proj= sqrt(x_proj*x_proj +y_proj*y_proj);

	{ //bound -pi< thet< pi
		tht_dum= fmod(tht_dum,2*M_PI);
		if(tht_dum<-M_PI) tht_dum+= +2*M_PI;
		if(tht_dum>M_PI) tht_dum+= -2*M_PI;
	}

	{pol_old[k]= pol[k]; tht_old[k]= tht[k]; zet_old[k]= zet_dum;}

	//record data
	rcrd_poin.push_back(k);
	rcrd_poin.push_back(rad_dum);
	rcrd_poin.push_back(tht_dum);
	rcrd_poin.push_back(r_proj);
	rcrd_poin.push_back(z_proj);
	const int numdat= 5; //number of data columns
	try{
		if(rcrd_poin.size()/numdat>MAX_POINPNTS){
			throw std::string("Max poincaré points reached: ")
				+std::to_string(MAX_POINPNTS);
		}
	}
	catch(const std::string& ex){
		std::cerr << ex <<'\n';
		return;
	}
}


void write_poin(){ // write poincare
	std::fstream poin;
	std::ofstream(outdir+"/poin.plt");
	poin.open(outdir+"/poin.plt",std::ios::app);
	const int numdat= 5;
	poin <<"#poin.plt    rad tht rcy zcy\n";
	for(int j=0;j<rcrd_poin.size()/numdat;j++){
		poin<<std::fixed << std::setw(4) <<std::setprecision(0)
			<< rcrd_poin[numdat*j+0] <<"  ";
		for(int i=1;i<numdat;i++){
			poin <<std::scientific <<std::setprecision(16)
				 << rcrd_poin[numdat*j+i] <<"  ";
			  }
		poin <<'\n';
	}
	poin.close();
	rcrd_poin.clear();
}
