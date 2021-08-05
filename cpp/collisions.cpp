//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#include "shared.h"
const double signnum(const double&);


const double dens_prof(
	const double& xd){ // density parabolic profile
	return 1.0- xd*xd;
}

const double temp_prof(
	const double& xd){ // temperature gaussian profile
	return exp(-xd);
}

const double deflection_factor(
	const double& x){ // deflection frequency factor
	return (1.0 -0.5/x/x)*erf(x)/x/x/x +exp(-x*x)/sqrt(M_PI)/x/x/x/x;
}

const double energy_factor(
	const double& x){ // energy frequency factor
	return 0.5*erf(x)/x/x/x -exp(-x*x)/sqrt(M_PI)/x/x;
}

const std::array<double,3> white(
	const double& sqx){ // White energy factor
	const double a_0= 0.75225;
	const double a_1= -0.238;
	const double a_2= 0.311;
	const double a_3= -0.0956;
	const double a_4= 0.0156;
	const double x= sqx*sqx;
	const double coeffs=
		a_0 +a_1*x +a_2*x*x +a_3*x*x*x +a_4*x*x*x*x;
	const double pfun=
		sqx*sqx*sqx*coeffs;
	const double psi_fun_approx=
		1.0 -1.0/(1.0 +pfun);
	const double psi_fun_prime=
		2*sqx*exp(-x)/sqrt(M_PI);
	const double psi_fun_prime_prime=
		exp(-x)/sqx/sqrt(M_PI)
		-2*sqx*exp(-x)/sqrt(M_PI);

	return {psi_fun_approx, psi_fun_prime, psi_fun_prime_prime};
}

void scatter(
	const int& k,
	const double& ran_num_1,
	const double& ran_num_2){

	const double ap0= 0.75225;
	const double ap1= -0.238;
	const double ap2= 0.311;
	const double ap3= -0.0956;
	const double ap4= 0.0156;

	const double prot= // test mass
		(double)test.i_mass;
	const double zprt= // test charge
		(double)test.i_chrg;
	const double massb= // background mass
		(double)plasma.i_mass;
	const double chgb= // background charge
		(double)plasma.i_chrg;

	// parameters
	const double ekmin= // minimum test kinetic energy in code units
		1.0e-2*KEV_TO_ERG*code.d_en/code.d_enax;
	double ekin= // test kinetic energy in code units
		en[k] -pt[k];
	const double eekev= // test kinetic energy in keV
		ekin*code.d_enax/KEV_TO_ERG;
	double dnb; // background density in cm^-3
	switch (plasma.dens_prof) {
		case DEN_FLAT: {
				dnb= plasma.d_dens;
				break;
			}
		case DEN_PARA: {
				dnb= plasma.d_dens
					*dens_prof(pol[k]/code.d_polw);
				break;
			}
	}
	const double dne= // electrons density
		dnb*chgb;
	double tikev; // background temperature in keV
	switch (plasma.temp_prof) {
		case TEM_FLAT: {
				tikev= (double)plasma.i_temp/1000;
				break;
			}
		case TEM_GAUS: {
				tikev= (double)plasma.i_temp/1000
					*temp_prof(pol[k]/code.d_polw);
				break;
			}
	}
	const double ti=
		tikev*KEV_TO_ERG/code.d_enax;
	const double tekev=
		tikev;
	const double te=
		tekev*KEV_TO_ERG/code.d_enax;
	const double vprt= // test particle velocity
		sqrt(2*ekin);
	const double vtb= // background velocity
		sqrt(2*ti/massb);
	const double vte= // electrons velocity
		sqrt(2*te*1836);

	const double cnst= // Coulomb logarithm units
		2.41e+11*(zprt/prot)*(zprt/prot)
		/(code.d_leng*code.d_freq*vprt)
		/(code.d_leng*code.d_freq*vprt)
		/(code.d_leng*code.d_freq*vprt);

	{ // background
		const double dnu0= // background species s^-1
			cnst*dnb*chgb*chgb;
		const double dumb= // velocity ratio
			vprt/vtb;
		const double dd= // dd= x, dumb= sqrt(x)
			dumb*dumb;
		const double d3=
			dumb*dumb*dumb;
		const double dum=
			d3*(ap0 +ap1*dd +ap2*dd*dd +ap3*dd*dd*dd +ap4*dd*dd*dd*dd);
		const double psib= // approx to psi fnctn
			1.0 -1.0/(1.0 +dum);
		const double fb=
			(2.0 -1.0/dd)*psib +2.257*dumb*exp(-dd);
		const double gb=
			2*chgb*psib/massb -2.257*dumb*exp(-dd);
		const double gbp=
			2.257*exp(-dd)*((1.0 +chgb/massb)*dumb -0.5/dumb);

		// find coll freq and convert from sec-1 to code units
		const double bmx=
			7.4e-3*sqrt(tikev*1.0e+13/(zprt*zprt*dnb));
		const double bmin1=
			1.4e-10*zprt*chgb/(eekev +tikev);
		const double massab=
			prot*massb/(prot +massb);
		const double bmin2=
			6.2e-10/(sqrt(massab*1836*tikev));
		const double bdum=
			std::max(bmin1,bmin2);
		const double clogb= // ib Coulomb log
			log(bmx/bdum);
		const double colb= // scatter
			abs(clogb*dnu0*fb)/code.d_freq;
		const double colbs= // slow
			abs(clogb*dnu0*gb)/code.d_freq;

		// scatter and drag
		{ // pitch angle
			const double dum=
				std::max(0.0,1.0 -ptch[k]*ptch[k]);
			const double psig=
				signnum(ran_num_1 -0.5)*sqrt(dum*colb*dt[k]*0.5);
			ptch[k]=
				ptch[k]*(1.0 -colb*dt[k]*0.5) +psig;
		}
		{ // energy
			const double esig=
				signnum(ran_num_2 -0.5)*sqrt(2*ekin*ti*colbs*dt[k]);
			const double efac0=
				ekin*(1.0 -massb*gbp/(gb*prot));
			ekin=
				ekin -colbs*dt[k]*efac0 +esig;
			ekin=
				std::max(ekin,ekmin); // truncate slowing down;
		}
	}

	{ // electrons
		const double enu0= // electrons  sec-1
			cnst*dne;
		const double dume=
			vprt/vte;
		const double dd=
			dume*dume;
		const double d3=
			dume*dume*dume;
		const double dum=
			d3*(ap0 +ap1*dd +ap2*dd*dd +ap3*dd*dd*dd +ap4*dd*dd*dd*dd);
		const double psie= // approx to psi fnctn
			1.0 -1.0/(1.0 +dum);
		const double fe=
			(2.0 -1.0/dd)*psie +2.257*dume*exp(-dd);
		const double ge=
			2*prot*psie*1836 -2.257*dume*exp(-dd);
		const double gep=
			2.257*exp(-dd)*((1.0 +prot*1836)*dume -0.5/dume);

		// find coll freq and convert from sec-1 to code units
		const double bmx=
			7.4e-3*sqrt(tekev*1.0e+13/(zprt*zprt*dne));
		const double bmin1=
			1.4e-10*zprt/(eekev +tekev);
		const double massab=
			prot/(1836*prot +1.0);
		const double bmin2=
			6.2e-10/(sqrt(massab*1836*tekev));
		const double bdum=
			std::max(bmin1,bmin2);
		const double cloge= // ie Coulomb log
			log(bmx/bdum);
		const double cole= // scatter
			abs(cloge*enu0*fe)/code.d_freq;
		const double coles= // slow
			abs(cloge*enu0*ge)/code.d_freq;

		// scatter and drag
		{ // pitch angle
			const double dum=
				std::max(0.0,1.0 -ptch[k]*ptch[k]);
			const double psig=
				signnum(ran_num_1 -0.5)*sqrt(dum*cole*dt[k]*0.5);
			ptch[k]=
				ptch[k]*(1.0 -cole*dt[k]*0.5) +psig;
		}
		{ // energy
			const double esig=
				signnum(ran_num_2 -0.5)*sqrt(2*ekin*te*coles*dt[k]);
			const double efac0=
				ekin*(1.0 -gep/(ge*1836*prot));
			ekin=
				ekin -coles*dt[k]*efac0 +esig;
			ekin=
				std::max(ekin,ekmin); // truncate slowing down;
		}
	}

	ptch[k]= std::min(1.0,ptch[k]);
	ptch[k]= std::max(-1.0,ptch[k]);
	en[k]= ekin +pt[k];

	rho[k]= ptch[k]*sqrt(2*ekin)/bf[k];
	mu[k]= ekin/bf[k] -0.5*rho[k]*rho[k]*bf[k];
}


void scatter_old(
	const int& k,
	const double& ran_num_1,
	const double& ran_num_2){
	// test particle parameters
	double mass_tst; // test mass in grams
	double chrg_tst; // test charge in statC
	switch (test.sign) {
		case CHG_IONS:
		{
			mass_tst= MASS_PROT*test.i_mass;
			chrg_tst= CHRG_PROT*test.i_chrg;
			break;
		}
		case CHG_ELEC:
		{
			mass_tst= MASS_ELEC;
			chrg_tst= CHRG_PROT;
			break;
		}
	}
	const double vel_tst=  // test particle speed in cm*s^-1
		sqrt(2*code.d_enax*en[k]/mass_tst);
	double ekin= en[k] -pt[k]; // test particle kinetic energy

	const int rand_sign_1= // random sign
		signnum(ran_num_1 -0.5);
	const int rand_sign_2= // random sign
		signnum(ran_num_2 -0.5);
	double pch_slow_sum= 0.0; // pitch angle slowdown term sum
	double pch_scat_sum= 0.0; // pitch angle scatter term sum
	double eng_slow_sum= 0.0; // energy slowdown term sum
	double eng_scat_sum= 0.0; // energy scatter term sum
for(int i=0;i<total_species;i++){
	switch (species[i].sign){
		case CHG_IONS:
		{
			// ions parameters
			const double mass_ion= // mass in g
				MASS_PROT*species[i].i_mass;
			const double chrg_ion= // charge in statC
				CHRG_PROT*species[i].i_chrg;
			double temp_ion; // temperature in ergs
			switch (species[i].temp_prof) {
				case TEM_FLAT: {
						temp_ion= KEV_TO_ERG*species[i].i_temp;
						break;
					}
				case TEM_GAUS: {
						temp_ion= KEV_TO_ERG*species[i].i_temp
							*temp_prof(pol[k]/code.d_polw);
						break;
					}
			}
			double dens_ion; // density in cm^-3
			switch (species[i].dens_prof) {
				case DEN_FLAT: {
						dens_ion= species[i].d_dens;
						break;
					}
				case DEN_PARA: {
						dens_ion= species[i].d_dens
							*dens_prof(pol[k]/code.d_polw);
						break;
					}
			}
			const double vel_ion=  // speed in cm*s^-1
				sqrt(3*temp_ion/mass_ion);
			const double x_ion= // velocity ratio
				vel_tst/vel_ion;

			const double bmin_ion= // impact parameter in cm
				chrg_tst*chrg_ion
				*(mass_tst +mass_ion)/mass_tst/mass_ion
				/vel_ion/vel_ion;
			const double bmax_ion= // debye length in cm
				sqrt(temp_ion/dens_ion/chrg_ion/chrg_ion);
			const double clog_ion= // coulomb logarithm
				log(bmax_ion/bmin_ion);
			const double nu_bi=    // braginskii frequency
				1.3333*sqrt(2*2*2*M_PI)
				*chrg_tst*chrg_tst*chrg_ion*chrg_ion/mass_tst/mass_ion
				*dens_ion*clog_ion
				/vel_ion/vel_ion/vel_ion/code.d_freq;

			// deflection operator
			const double coef_di= // deflection coefficient
				1.5*sqrt(0.5*M_PI)
				*chrg_tst*chrg_tst*mass_ion*mass_ion
				/chrg_ion/chrg_ion/mass_tst/mass_tst;
			const double nu_di= // deflection frequency
					coef_di*nu_bi*deflection_factor(x_ion);
			pch_slow_sum+= // pitch angle slowdown term
				-nu_di*ptch[k]*dt[k];
			pch_scat_sum+= // pitch angle scatter term
				(1.0- ptch[k]*ptch[k])*nu_di*dt[k];

			// energy operator
			const double coef_ei= // energy coefficient
				3*sqrt(M_PI/2)
				*mass_ion*chrg_tst*chrg_tst
				/chrg_ion/chrg_ion/mass_tst;
			const double nu_ei= // energy frequency
				coef_ei*nu_bi*energy_factor(x_ion);
			// const double dum=
			// 	ekin -1.5*mass_ion*exp(-x_ion*x_ion)*temp_ion
			// 	/energy_factor(x_ion)/mass_tst/sqrt(M_PI)/code.d_enax;
			const std::array<double,3> psi_fun= white(x_ion);
			const double dum=
				mass_ion*(psi_fun[0] +psi_fun[1])/mass_tst
				-psi_fun[1] -psi_fun[2];
			eng_slow_sum+= // energy slowdown term
				-2*nu_ei*ekin*dum*dt[k];
			eng_scat_sum+= // energy scatter term
				nu_ei*temp_ion*ekin*dt[k]/code.d_enax;
			break;
		}
		case CHG_ELEC:
		{
			// electrons parameters
			const double mass_ele= // mass in g
				MASS_ELEC*species[i].i_mass;
			const double chrg_ele= // charge in statC
				CHRG_PROT*abs(species[i].i_chrg);
			double temp_ele; // temperature in ergs
			switch (species[i].temp_prof) {
				case TEM_FLAT: {
						temp_ele= KEV_TO_ERG*species[i].i_temp;
						break;
					}
				case TEM_GAUS: {
						temp_ele= KEV_TO_ERG*species[i].i_temp
							*temp_prof(pol[k]/code.d_polw);
						break;
					}
			}
			double dens_ele; // density in cm^-3
			switch (species[i].dens_prof) {
				case DEN_FLAT: {
						dens_ele= species[i].d_dens;
						break;
					}
				case DEN_PARA: {
						dens_ele= species[i].d_dens
							*dens_prof(pol[k]/code.d_polw);
						break;
					}
			}
			const double vel_ele=  // speed in cm*s^-1
				sqrt(3*temp_ele/mass_ele);
			const double x_ele= // velocity ratio
				vel_tst/vel_ele;

			const double bmin_ele= // electrons impact parameter in cm
				chrg_tst*chrg_ele
				*(mass_tst +mass_ele)/mass_tst/mass_ele
				/vel_ele/vel_ele;
			const double bmax_ele= // electrons debye length in cm
				sqrt(temp_ele/dens_ele/chrg_ele/chrg_ele);
			const double clog_ele= // electrons coulomb logarithm
				log(bmax_ele/bmin_ele);
			const double nu_be= // braginskii electrons frequency
				5.3333*sqrt(M_PI)
				*chrg_tst*chrg_tst*chrg_ele*chrg_ele/mass_ele/mass_ele
				*dens_ele*clog_ele
				/vel_ele/vel_ele/vel_ele/code.d_freq;

			// deflection operator
			const double coef_de=
				0.75*sqrt(M_PI)
				*chrg_tst*chrg_tst*mass_ele*mass_ele
				/chrg_ele/chrg_ele/mass_tst/mass_tst;
			const double nu_de=
				coef_de*nu_be*deflection_factor(x_ele);
			pch_slow_sum+= // pitch angle slowdown term
				-nu_de*ptch[k]*dt[k];
			pch_scat_sum+= // pitch angle scatter term
				(1.0- ptch[k]*ptch[k])*nu_de*dt[k];

			// energy operator
			const double coef_ee=
				1.5*sqrt(M_PI)
				*mass_ele*chrg_tst*chrg_tst/mass_tst/chrg_ele/chrg_ele;
			const double nu_ee=
				coef_ee*nu_be*energy_factor(x_ele);
			// const double dum=
			// 	ekin -1.5*mass_ele*exp(-x_ele*x_ele)*temp_ele
			// 	/energy_factor(x_ele)/mass_tst/sqrt(M_PI)/code.d_enax;
			const std::array<double,3> psi_fun= white(x_ele);
			const double dum=
				mass_ele*(psi_fun[0] +psi_fun[1])/mass_tst
				-psi_fun[1] -psi_fun[2];
			eng_slow_sum+= // energy slowdown term
				-2*nu_ee*ekin*dum*dt[k];
			eng_scat_sum+= // energy scatter term
				nu_ee*temp_ele*ekin/code.d_enax*dt[k];
			break;
		}
	}
}

	{ // pitch-angle collision step
		pch_scat_sum= // scatter term
			rand_sign_1*sqrt(pch_scat_sum);
		ptch[k]+= pch_slow_sum +pch_scat_sum;
		ptch[k]= std::min(1.0,ptch[k]);
		ptch[k]= std::max(-1.0,ptch[k]);
	}

	{ // energy collision step
		const double ekin_min= // minimum energy bound
			1.0e-2*1*KEV_TO_ERG/code.d_enax;
		eng_scat_sum= // scatter term
			rand_sign_2*2*sqrt(eng_scat_sum);
		ekin+= eng_slow_sum +eng_scat_sum;
		ekin= std::max(ekin,ekin_min); // truncate slowdown
		en[k]= ekin;
	}

	{ // update variables
		rho[k]= ptch[k]*sqrt(2*en[k])/bf[k];
		mu[k]= en[k]/bf[k] -0.5*rho[k]*rho[k]*bf[k];
	}
}
