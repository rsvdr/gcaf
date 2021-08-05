//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#include "shared.h"
#include "define.h"
#include <chrono>
#include <omp.h>
#include <time.h>

// Variable declarations
std::array<double,MAX_SDEL*MAX_TRAN> rand_list_1;
std::array<double,MAX_SDEL*MAX_TRAN> rand_list_2;

// Function declarations
void info();
void boost_parse(int ac, char* av[]);
void check_bounds();
void inivar();
void set(const int&);
void logger();
void record_traj(const int&, const int&);
void record_flux(const int&, const int&);
void record_mome(const int&, const int&);
void record_dist(const int&, const int&);
void record_poin(const int&);
void write_traj();
void write_flux();
void write_mome();
void write_dist();
void write_poin();
void onestep(const int&);
void perturbation(const int&);
void scatter(const int&, const double&, const double&);
void update(const int&);
const double urand();


struct Timer{ // cpu timer
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<float> duration;
	Timer(){
		start= std::chrono::high_resolution_clock::now();
	}
	~Timer(){
		end= std::chrono::high_resolution_clock::now();
		duration= end -start;
		float dur = duration.count();
		std::cout <<"Timer: " <<dur <<" s"<<'\n';
	}
};

void traj_loop(){ // trajectory simulation loop
	srand(time(0));
	int step_time= 0; // time step counter
	int skip_traj= 0; // traj record step counter
	for(int i=0;i<sim.i_tran*sim.i_sdel;i++){ // random numbers
		rand_list_1[i]= (double)rand()/RAND_MAX;
		rand_list_2[i]= (double)rand()/RAND_MAX;
	}
	const int k= 0;
	for(int i=0;i<sim.i_sdel*sim.i_tran;i++){ // time loop
		onestep(k);
		update(k);
		if(swt.coll)
			scatter(k,rand_list_1[step_time],rand_list_2[step_time]);
		if(swt.pert) perturbation(k);
		if(tim[k]>skip_traj*code.d_trun/sim.i_traj)
			{record_traj(k,skip_traj); skip_traj++;}
		step_time++;
	}
}

void ense_loop(){ // ensemble simulation loop
	srand(time(0));
	// #pragma omp parallel for
	for(int k=0;k<sim.i_part;k++){ //particle loop
		int thread_id = omp_get_thread_num();
		int step_time= 0; // time step counter
		int skip_mome= 0; // mome record step counter
		int skip_dist= 0; // dist record step counter
		flux[k]= 0;
		for(int i=0;i<sim.i_tran*sim.i_sdel;i++){ // random numbers
			rand_list_1[i]= (double)rand()/RAND_MAX;
			rand_list_2[i]= (double)rand()/RAND_MAX;
		}
		// for(int i=0;i<sim.i_sdel*sim.i_tran/sim.i_dist;i++)
		// { // initial time loop
		// 	if(swt.pert) perturbation(k);
		// 	onestep(k);
		// 	update(k);
		// }
		record_dist(k,skip_dist); // initial distribution
		skip_dist++;
		for(int i=0;i<sim.i_sdel*sim.i_tran;i++){ // main time loop
			if(swt.pert) perturbation(k);
			if(swt.coll)
				scatter(k,rand_list_1[step_time],rand_list_2[step_time]);
			onestep(k);
			if(swt.flux) record_flux(k,skip_mome);
			if(tim[k]>skip_mome*code.d_trun/sim.i_mome)
				{record_mome(k,skip_mome); skip_mome++;}
			if(tim[k]>skip_dist*code.d_trun/sim.i_dist)
				{record_dist(k,skip_dist); skip_dist++;}
			update(k);
			step_time++;
		}
		record_dist(k,skip_dist); // final distribution
		if(swt.verb)
			printf("Particle %d finished in thread %d\n",k,thread_id);
	}
}

void poin_loop(){ // poincare simulation loop
	// #pragma omp parallel for
	for(int k=0;k<sim.i_part;k++){ // particle loop
		for(int i=0;i<sim.i_sdel*sim.i_tran;i++){ // time loop
			if(swt.pert) perturbation(k);
			onestep(k);
			update(k);
			#pragma omp critical
			{record_poin(k);}
		}
		if(swt.verb) printf("Particle %d finished\n",k);
	}
}

int main(int argc,char* argv[]){
	Timer timer;
try{
	info();
	boost_parse(argc,argv);
	check_bounds();
	inivar();
	for(int k=0;k<sim.i_part;k++) set(k);
	std::clog <<"Variables initialized successfully...\n";
	logger();

#if 1
	switch (swt.simu){
		case SIM_TRAJ: {
			std::clog <<"Trajectory simulation start...\n";
			traj_loop();
			write_traj();
			std::clog <<"Finished successfully!\n";
			break;
		}
		case SIM_ENSE: {
			std::clog <<"Ensemble simulation start...\n";
			ense_loop();
			if(swt.flux) write_flux();
			write_mome();
			write_dist();
			std::clog <<"Finished successfully!\n";
			break;
		}
		case SIM_POIN: {
			std::clog <<"Poincaré simulation start...\n";
			poin_loop();
			write_poin();
			std::clog <<"Finished successfully!\n";
			break;
		}
	}
#endif
}
catch(const std::string& ex){
	std::cerr <<"[ERROR] " <<ex <<'\n';
	return 0;
}
catch(const int& ex){
	return 0;
}
catch(std::exception& ex){
	std::cerr <<"[ERROR] " <<ex.what() <<'\n';
	return 0;
}
}
