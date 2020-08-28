//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://github.com/Saavestro/gcafpp //
//----------------------------------------------------------------//
#include "shared.h"


const double normal_dist(const double&, const double&);
const double signnum(const double&);
void field_tkmk(const int&);
void field_stel(const int&);
void electrostatic_potential(const int&);
void q_currents(const int&);
const double urand();


void inivar(){ //initialize simulation variables
	if(swt.simu==SIM_TRAJ) sim.i_part= 1;
	code.d_leng= // characteristic length
		field.i_rmin;
	code.d_poli= // initial poloidal flux
		0.5*pow(test.d_radi*field.i_rmin/code.d_leng,2);
	code.d_polw= // on wall poloidal flux
		0.5*pow(1.0*field.i_rmin/code.d_leng,2);
	code.d_thet= // initial theta
		test.i_thet*M_PI/180;
	code.d_zeta= // initial theta
		test.i_zeta*M_PI/360;
	code.d_ptch= // initial pitch angle
		cos(test.i_ptch*M_PI/180);

	if(test.sign==CHG_IONS){
		code.d_freq= // characteristic frequency
			CHRG_PROT*TSLA_TO_GAUS*test.i_chrg*field.i_baxi
			/test.i_mass/MASS_PROT/CSPEED;
		code.d_enax= // characteristic'energy
			MASS_PROT*test.i_mass*pow(code.d_leng*code.d_freq,2);
	}
	else{
		code.d_freq= // characteristic frequency
		 	1836.1526*CHRG_PROT*TSLA_TO_GAUS*test.i_chrg*field.i_baxi
			/test.i_mass/MASS_PROT/CSPEED;
		code.d_enax= // characteristic energy
			MASS_PROT*test.i_mass*pow(code.d_leng*code.d_freq,2)/1836.1526;
	}
	code.d_pert_freq= // perturbation frequency
		field.i_freq/code.d_freq;
	code.d_pt= // electrostatic potential
		KEV_TO_ERG*field.d_epot/code.d_enax;

	const int bt= field.i_baxi;
	code.d_en= test.d_ener;
	const double go= field.i_rmaj*bt/code.d_leng/field.i_baxi;
	const double tran= 2*M_PI*go/bt/sqrt(2*KEV_TO_ERG*code.d_en/code.d_enax);
	code.d_trun= sim.i_tran*tran;
	code.d_dt= tran/sim.i_sdel;

	code.d_zeta_poin= sim.i_surf*2*M_PI/360 +2*M_PI/1e4;
	{ //bound 0 < code.d_zeta_poin < 2*pi
		code.d_zeta_poin= fmod(code.d_zeta_poin,2*M_PI);
		if(code.d_zeta_poin<0.0) code.d_zeta_poin+= +2*M_PI;
	}

	// collision parameters
	switch (test.sign) {
		case CHG_IONS:
		{
			code.d_test_mass= MASS_PROT*test.i_mass;
			code.d_test_chrg= CHRG_PROT*test.i_chrg;
			break;
		}
		case CHG_ELEC:
		{
			code.d_test_mass= MASS_ELEC;
			code.d_test_chrg= CHRG_PROT;
			break;
		}
	}
	plasma.i_mass= 1;
	plasma.i_chrg= 1;
}

double set_pitch(pch_dist_t str){
	switch(str){
		case PCH_PAR:
			return 1.0;
			break;
		case PCH_ANT:
			return -1.0;
			break;
		case PCH_ISO:
			return signnum(0.5 -urand())*urand();
			break;
		case PCH_FIX:
			return code.d_ptch;
			break;
	}
}


void set(const int& k){ //set test particle variables
	//set pol,thet,zet,en,ptch
	switch(swt.simu){
		case SIM_TRAJ: {
				pol[k]= code.d_poli;
				tht[k]= code.d_thet;
				zet[k]= code.d_zeta;
				ptch[k]= set_pitch(test.ptch_dist);
				en[k]= KEV_TO_ERG*code.d_en/code.d_enax;
				break;
				}
		case SIM_ENSE: {
				pol[k]= code.d_poli;
				tht[k]= code.d_thet; //signnum(0.5 -urand())*M_PI*urand();
				zet[k]= code.d_zeta; //2*M_PI*urand();
				ptch[k]= set_pitch(test.ptch_dist);
				double endum= code.d_en;
				{ // maxwellian distribution
					const double vx= sqrt(KEV_TO_ERG*endum/MASS_PROT/test.i_mass)*normal_dist(0.0,1.0);
					const double vy= sqrt(KEV_TO_ERG*endum/MASS_PROT/test.i_mass)*normal_dist(0.0,1.0);
					const double vz= sqrt(KEV_TO_ERG*endum/MASS_PROT/test.i_mass)*normal_dist(0.0,1.0);
					const double vsq= vx*vx +vy*vy +vz*vz;
					if(test.ener_dist==ENG_MAXW)
						endum= MASS_PROT*test.i_mass*vsq/2/KEV_TO_ERG;
				}
				en[k]= KEV_TO_ERG*endum/code.d_enax;
				break;
				}
		case SIM_POIN: {
				const double rone_dum= sim.d_rlow*field.i_rmin;
				const double rtwo_dum= sim.d_rupp*field.i_rmin;
				const double rdum= rone_dum+ k*(rtwo_dum- rone_dum)/sim.i_part;
				pol[k]= 0.5*rdum*rdum/code.d_leng/code.d_leng;
				tht[k]= code.d_thet;
				zet[k]= code.d_zeta;
				ptch[k]= set_pitch(test.ptch_dist);
				en[k]= KEV_TO_ERG*code.d_en/code.d_enax;
				break;
				}
	}

	//set t,dt,rho,mu,pz
	q_currents(k);
	if(swt.conf==DEV_TKMK){
		field_tkmk(k); //needs pol,tht,zet already set
		electrostatic_potential(k);
	}
	else{
		field_stel(k);
		electrostatic_potential(k);
	}
	tim[k]= 0.0;
	dt[k]= code.d_dt;
	rho[k]= ptch[k]*sqrt(2.0*en[k])/bf[k];
	mu[k]= en[k]/bf[k]- .5*rho[k]*rho[k]*bf[k];
	pz[k]= cg[k]*(rho[k] +al[k]) -pol[k];
}
