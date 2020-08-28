//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://github.com/Saavestro/gcafpp //
//----------------------------------------------------------------//
#include "shared.h"

const double temp_prof(const double&);


void q_currents(const int &k){ // q profile and currents
	const double ew= // Inverse aspect ratio
		(double)field.i_rmin/field.i_rmaj;
	const double ra= // r position over minor radius
		sqrt(pol[k]/code.d_polw);
	const double rap= // ra pol derivative
		0.5*ra/pol[k];
	const double lmsq= // lambda squared
		(double)field.i_rmin*field.i_rmin*field.i_qaxi
		/(field.i_qwll -field.i_qaxi);

	qf[k]= // q profile
		(double)field.i_qaxi
		+field.i_qaxi*ra*ra*field.i_rmin*field.i_rmin/lmsq;
	qf_pol[k]= // q profile derivative
		field.i_qaxi*field.i_rmin*field.i_rmin/lmsq/code.d_polw;
	cg[k]= // poloidal current
		field.i_rmaj/(1.0 +ra*ew*cos(tht[k]))/code.d_leng;
	cg_pol[k]= // poloidal current derivative
		rap*ew*cos(tht[k])/(1.0 +ra*ew*cos(tht[k]))/code.d_leng;
	ci[k]= // toroidal current
		field.i_rmin*ra*ew/(1.0 +ra*ew*cos(tht[k]))/qf[k]/code.d_leng;
	const double ci_polft= rap*ew;
	const double ci_polst= -qf_pol[k]/qf[k]/qf[k];
	const double ci_poltt= rap*ew*cos(tht[k])/(1.0 +ra*ew*cos(tht[k]));
	ci_pol[k]= // toroidal current derivative
		ci_polft*ci_polst*ci_poltt/code.d_leng;
}


void field_tkmk(const int& k){ // tokamak magnetic field
	q_currents(k);
	const double ew= (double)field.i_rmin/field.i_rmaj;
	const double ra= sqrt(pol[k]/code.d_polw);
	const double rap= 0.5*ra/pol[k];
	const double bft= (double)field.i_baxi;
	const double bst= ra*ew*field.i_baxi/qf[k];
	const double bstp= rap*ew*field.i_baxi/qf[k]
		-ra*ew*field.i_baxi*qf_pol[k]/qf[k]/qf[k];

	bf[k]= // b field
		sqrt(bft*bft +bst*bst)/field.i_baxi/(1.0 +ew*ra*cos(tht[k]));
	const double bfpft= - ew*rap*cos(tht[k])/(1.0 +ew*ra*cos(tht[k]))
		/(1.0 +ew*ra*cos(tht[k]));
	const double bfpst= bst*bstp/sqrt(bft*bft +bst*bst);
	bf_pol[k]= // b field pol derivative
		bfpst/field.i_baxi/(1.0 +ew*ra*cos(tht[k]))
		+bfpft*sqrt(bft*bft +bst*bst)/field.i_baxi;
	bf_tht[k]= // b field tht derivative
		ew*ra*sin(tht[k])*sqrt(bft*bft +bst*bst)
		/field.i_baxi/(1.0 +ew*ra*cos(tht[k]))/(1.0 +ew*ra*cos(tht[k]));
	bf_zet[k]= // b field zet derivative
		0.0;
}


void electrostatic_potential(const int& k){ // electrostatic potential
	if(field.potp==POT_CONS){
		const double dum= pol[k]/code.d_polw;
		pt[k]= code.d_pt*pol[k]*pol[k]/code.d_polw/code.d_polw;
		pt_pol[k]= code.d_pt*pol[k]/code.d_polw/code.d_polw;
		pt_tht[k]= 0.0;
		pt_zet[k]= 0.0;
	}
	else{ // needs to be fixed
		const double tmpdum= KEV_TO_ERG*species[0].i_temp/code.d_enax;
		pt[k]= code.d_pt +tmpdum*temp_prof(pol[k]/code.d_polw);
		pt_pol[k]= code.d_pt
			-tmpdum*temp_prof(pol[k]/code.d_polw)/code.d_polw;
		pt_tht[k]= 0.0;
		pt_zet[k]= 0.0;
	}
}


void field_stel(const int& k){ // stellarator magnetic field
	q_currents(k);
	const int lc= 2;
	const int mc= 6;
	const double sg= -1;
	const double ew0= (double)field.i_rmin/field.i_rmaj;
	const double ew1= 0.4;
	const double ewm= 0.4;
	const double ewp= 0.4;
	const double pfrs0= sqrt(pol[k]/code.d_polw);
	const double pfrs1= pow(pfrs0,lc);
	const double pfrsm= pow(pfrs0,lc-1);
	const double pfrsp= pow(pfrs0,lc+1);
	const double pfrs0p= 0.5*pfrs0/pol[k];
	const double pfrs1p= 0.5*lc*pfrs1/pol[k];
	const double pfrsmp= 0.5*(lc -1)*pfrsm/pol[k];
	const double pfrspp= 0.5*(lc +1)*pfrsp/pol[k];
	const double coss= cos(tht[k]);
	const double sins= sin(tht[k]);
	const double cosl= cos(lc*tht[k] -mc*zet[k]);
	const double sinl= sin(lc*tht[k] -mc*zet[k]);
	const double coslm= cos((lc -1)*tht[k] -mc*zet[k]);
	const double sinlm= sin((lc -1)*tht[k] -mc*zet[k]);
	const double coslp= cos((lc +1)*tht[k] -mc*zet[k]);
	const double sinlp= sin((lc +1)*tht[k] -mc*zet[k]);

	bf[k]= 1.0 -ew0*pfrs0*coss -ew1*pfrs1*cosl
		+0.5*sg*(ewm*pfrsm*coslm +ewp*pfrsp*coslp);
	bf_pol[k]= -ew0*pfrs0p*coss -ew1*pfrs1p*cosl
		+0.5*sg*(ewm*pfrsmp*coslm +ewp*pfrspp*coslp);
	bf_tht[k]= ew0*pfrs0*sins +lc*ew1*pfrs1*sinl
		-0.5*sg*((lc -1)*ewm*pfrsm*sinlm +(lc +1)*ewp*pfrsp*sinlp);
	bf_zet[k]= -mc*ew1*pfrs1*sinl
		+0.5*sg*(mc*ewm*pfrsm*sinlm +mc*ewp*pfrsp*sinlp);
}


void perturbation(const int& k){ // magnetic perturbation
	const double ph= (double)field.i_phas*M_PI/180;
	const int mi= field.i_mmod;
	const int mf= field.i_mmod;
	const int n= field.i_nmod;
	const int rp= n;
	for(int m=mi;m<=mf;m++){
		double facx= (double)m/(m +rp);
		double ra= field.d_ampx*pow(pol[k]/code.d_polw/facx,m)
			*pow(1.0 -pol[k]/code.d_polw,rp)/pow(1.0 -facx,rp);
		double rap= m*ra/pol[k] -rp*ra/code.d_polw/(1.0 -pol[k]/code.d_polw);
		double sin_long= sin(n*zet[k] -m*tht[k]
				+ph -code.d_pert_freq*tim[k]);
		double cos_long= cos(n*zet[k] -m*tht[k]
				+ph -code.d_pert_freq*tim[k]);
		al[k]= // alpha perturbation
			ra*sin_long;
		al_pol[k]= // alpha pol derivative
			rap*sin_long;
		al_tht[k]= // alpha tht derivative
			m*ra*cos_long;
		al_zet[k]= // alpha zet derivative
			n*ra*cos_long;
		al_tim[k]= // alpha tim derivative
			code.d_pert_freq*ra*cos_long;
	}
}
