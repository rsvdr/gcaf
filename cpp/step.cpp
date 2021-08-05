//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#include "shared.h"

void q_currents(const int&);
void field_tkmk(const int&);
void field_stel(const int&);
void electrostatic_potential(const int&);


void onestep(const int& k){ // integration step
	double eq[4]; // differential equation
	double f0[4]; // RK4 first step
	double f1[4]; // RK4 second step
	double f2[4]; // RK4 third step
	double f3[4]; // RK4 fourth step
	double u1[4]= {pol[k], tht[k], zet[k], rho[k]};
	double u0[4]= {u1[0], u1[1], u1[2], u1[3]};

	for(int i=0;i<4;i++){ //RK4 integrator
		q_currents(k);
		if(swt.conf==DEV_TKMK){
			field_tkmk(k);
			electrostatic_potential(k);
		}
		else{
			field_stel(k);
			electrostatic_potential(k);
		}
		const double rbb= u1[3]*bf[k]*bf[k];
		const double dedb= mu[k] +u1[3]*u1[3]*bf[k];
		const double det= cg[k]*qf[k] +ci[k]
			+(test.sign*u1[3] +al[k])*(cg[k]*ci_pol[k] -ci[k]*cg_pol[k]);
		const double fac1=
			1.0 -cg_pol[k]*(test.sign*u1[3] +al[k]) -cg[k]*al_pol[k];
		const double fac2=
			qf[k] +ci_pol[k]*(test.sign*u1[3] +al[k]) +ci[k]*al_pol[k];
		const double fac3=
			test.sign*ci[k]*al_zet[k]- test.sign*cg[k]*al_tht[k];
		eq[0]= // pol equation
			-test.sign*cg[k]*dedb*bf_tht[k]/det
			+test.sign*ci[k]*dedb*bf_zet[k]/det
			+cg[k]*rbb*al_tht[k]/det
			-ci[k]*rbb*al_zet[k]/det
			-test.sign*cg[k]*pt_tht[k]/det
			+test.sign*ci[k]*pt_zet[k]/det;
		eq[1]= // tht equation
			test.sign*cg[k]*dedb*bf_pol[k]/det
			+rbb*fac1/det
			+test.sign*cg[k]*pt_pol[k]/det;
		eq[2]= // zet equation
			-test.sign*ci[k]*dedb*bf_pol[k]/det
			+rbb*fac2/det
			-test.sign*ci[k]*pt_pol[k]/det;
		eq[3]= // rho equation
			-fac1*dedb*bf_tht[k]/det
			-fac2*dedb*bf_zet[k]/det
			+fac3*dedb*bf_pol[k]/det
			-fac1*pt_tht[k]/det
			-fac2*pt_zet[k]/det
			+fac3*pt_pol[k]/det
			-al_tim[k];

		switch(i){
		case 0: { // RK4 first step
				for(int j=0;j<4;j++)
					{f0[j]= eq[j]; u1[j]= u0[j] +0.5*dt[k]*f0[j];}
				break;
				}
		case 1: { // RK4 second step
				for(int j=0;j<4;j++)
					{f1[j]= eq[j]; u1[j]= u0[j] +0.5*dt[k]*f1[j];}
				break;
				}
		case 2: { // RK4 third step
				for(int j=0;j<4;j++)
					{f2[j]= eq[j]; u1[j]= u0[j] +0.5*dt[k]*eq[j];}
				break;
				}
		case 3: { // RK4 fourth step
				for(int j=0;j<4;j++)
					{f3[j]= eq[j];
					u1[j]= u0[j] +dt[k]*(f0[j] +2*f1[j] +2*f2[j] +f3[j])/6;}
				break;
				}
		}
	}
	pol[k]= u1[0];
	tht[k]= u1[1];
	zet[k]= u1[2];
	rho[k]= u1[3];
}

void update(const int& k){ // update particle variables
	q_currents(k);
	if(swt.conf==DEV_TKMK) field_tkmk(k);
	else field_stel(k);
	toten[k]= 0.5*bf[k]*bf[k]*rho[k]*rho[k] +mu[k]*bf[k] +pt[k];
	en[k]= toten[k] -pt[k];
	ptch[k]= rho[k]*bf[k]/sqrt(2*en[k]);
	tim[k]+= dt[k];
	pz[k]= cg[k]*(rho[k] +al[k]) -pol[k];
	rad_prev[k]= code.d_leng*sqrt(2*pol[k]);
}
