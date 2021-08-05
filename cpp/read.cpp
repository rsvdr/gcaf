//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://git.rsaavedra.xyz/?p=gcaf.git //
//----------------------------------------------------------------//
#include "shared.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

std::vector<std::string> strip(std::string);
po::variables_map vm;

void read_slice(){ // read slices data
	std::vector<double> t_dat;
	std::vector<double> r_dat;
	std::vector<std::string> rt_tokens;
	std::string data_line;
	std::fstream surf;
	for(int j=0;j<total_slices;j++){
		const std::string slice_file=
			vm["slices"].as<std::vector<std::string>>()[j];
		surf.open(slice_file,std::ios::in);
		while(std::getline(surf,data_line)){
			rt_tokens =strip(data_line);
			t_dat.push_back(std::stod(rt_tokens[0]));
			r_dat.push_back(std::stod(rt_tokens[1]));
		}

		ip_size= t_dat.size();
		for(int i=0;i<t_dat.size();i++)
			{tht_ip[j+ i]= t_dat[i]; rad_ip[j+ i]= r_dat[i];}
		t_dat.clear();
		r_dat.clear();
		surf.close();
	}
}


struct coordinate { // parse values as "first second"
	std::string coord_first;
	std::string coord_second;
};
void validate(boost::any& v,
	const std::vector<std::string>& values,
	coordinate*, int) {
	coordinate c;
	std::vector<std::string> dvalues;
	for(const auto& val : values){
		std::stringstream ss(val);
		std::copy(std::istream_iterator<std::string>(ss),
			std::istream_iterator<std::string>(),
			std::back_inserter(dvalues));
		if(!ss.eof())
			throw po::invalid_option_value("Invalid specification");
	}
	if(dvalues.size()<1){
		throw po::invalid_option_value("Invalid specification");
	}
	c.coord_first= dvalues[0];
	if(dvalues.size()==2) c.coord_second = dvalues[1];
	v= c;
}

void boost_parse(int ac, char* av[]){ // boost library parsing
	using namespace std;
	std::string indir;  // input file

	// Command line options
	po::options_description command_options("Command line options");
	command_options.add_options()
		("help,h",
			"help message")
		("input-file,i",po::value<string>(&indir)->required(),
			"input file")
		("output-path,o",po::value<string>(&outdir)->required(),
			"output path")
		("simulation_type,t",po::value<string>()->required(),
			"simulation type traj/ense/poin")
		("device,d",po::value<string>(),
			"magnetic field type tkmk/stel")
		("verbose,v",
			"verbose")
		("screen-log,l",
			"show log to screen")
		("collisions,c",
			"collisions")
		("perturbation,p",
			"perturbation")
		("density",po::value<double>(),
			"plasma density in cm^-3")
	    ("slices,s",po::value<vector<string>>()->multitoken(),
	        "slices data [path, number-of-files]")
	    ;
	// Simulation options
	po::options_description config_simulation("Simulation options");
	config_simulation.add_options()
	    ("simulation.number_of_particles",
			po::value<int>(&sim.i_part)->multitoken(),
			"number of particles")
		("simulation.toroidal_transits",
			po::value<int>(&sim.i_tran)->multitoken(),
			"number of toroidal transits")
		("simulation.step_delta",
			po::value<int>(&sim.i_sdel)->multitoken(),
			"one toroidal transit/N")
		("simulation.trajectory_points",
			po::value<int>(&sim.i_traj)->multitoken(),
			"number of points in traj.plt")
		("simulation.moment_points",
			po::value<int>(&sim.i_mome)->multitoken(),
			"number of points in mome.plt")
		("simulation.distribution_points",
			po::value<int>(&sim.i_dist)->multitoken(),
			"number of distribution snapshots")
		("simulation.poincare_bounds",
			po::value<coordinate>()->multitoken(),
			"[lower bound in r/rmin, lower bound in r/rmin]")
		("simulation.poincare_plane",
			po::value<int>(&sim.i_surf)->multitoken(),
			"constant zeta plane in grades")
	;
	// Field options
	po::options_description config_field("Field options");
	config_field.add_options()
		("field.field_magnitude",
			po::value<int>(&field.i_baxi)->multitoken(),
			"on axis magnetic field magnitude")
		("field.q_on_axis",
			po::value<int>(&field.i_qaxi)->multitoken(),
			"on axis q profile")
		("field.q_on_wall",
			po::value<int>(&field.i_qwll)->multitoken(),
			"on wall q profile")
		("field.major_radius",
			po::value<int>(&field.i_rmaj)->multitoken(),
			"major radius in cm")
		("field.minor_radius",
			po::value<int>(&field.i_rmin)->multitoken(),
			"minor radius in cm")
		("field.perturbation_amplitude",
			po::value<double>(&field.d_ampx)->multitoken(),
			"perturbation amplitude in")
		("field.perturbation_frequency",
			po::value<int>(&field.i_freq)->multitoken(),
			"perturbation frequency in Hz")
		("field.perturbation_phase",
			po::value<int>(&field.i_phas)->multitoken(),
			"perturbation phase in grades")
		("field.m_mode",
			po::value<int>(&field.i_mmod)->multitoken(),
			"toroidal mode number")
		("field.n_mode",
			po::value<int>(&field.i_nmod)->multitoken(),
			"poloidal mode number")
		("field.electrostatic_potential",
			po::value<coordinate>()->multitoken(),
			"[profile \'cons\' or \'pres\', electrostatic potential in keV]")
		;
	// Plasma options
	po::options_description config_plasma("Plasma options");
	config_plasma.add_options()
		("plasma.temperature",
			po::value<coordinate>()->multitoken(),
			"[profile \'flat\' or \'gaus\', temperature in eV]")
		("plasma.density",
			po::value<coordinate>()->multitoken(),
			"[profile \'flat\' or \'para\', density in cm^-3]")
		;
	// impurities options
	po::options_description config_impurities("Impurities options");
	config_impurities.add_options()
		("impurities.type",
			po::value<vector<string>>()->multitoken(),
			"particle type \'ions\' or \'electrons\'")
		("impurities.mass",
			po::value<vector<int>>()->multitoken(),
			"mass in protons")
		("impurities.charge",
			po::value<vector<int>>()->multitoken(),
			"charge in protons")
		("impurities.temperature",
			po::value<vector<coordinate>>()->multitoken(),
			"[profile \'flat\' or \'gaus\', temperature in keV]")
		("impurities.density",
			po::value<vector<coordinate>>()->multitoken(),
			"[profile \'flat\' or \'para\', density in cm^-3]")
		;
	// Test particle options
	po::options_description config_particle("Test particle options");
	config_particle.add_options()
		("particle.type_test",
			po::value<string>()->multitoken(),
			"particle type \'ions\' or \'electrons\'")
		("particle.mass_test",
			po::value<int>(&test.i_mass)->multitoken(),
			"test mass in protons")
		("particle.charge_test",
			po::value<int>(&test.i_chrg)->multitoken(),
			"test charge in protons")
		("particle.radial_position",
			po::value<double>(&test.d_radi)->multitoken(),
			"radial position in minor radius")
		("particle.theta",
			po::value<int>(&test.i_thet)->multitoken(),
			"theta coordinate in grades")
		("particle.zeta",
			po::value<int>(&test.i_zeta)->multitoken(),
			"zeta coordinate in grades")
		("particle.pitch_angle",
			po::value<coordinate>()->multitoken(),
			"[distribution \'par\', \'ant\', \'iso\' or \'fix\', pitch angle in grades]")
		("particle.energy",
			po::value<coordinate>()->multitoken(),
			"[distribution \'mono\' or \'maxw\', energy in keV]")
		;
	// Group options
	po::options_description global_options;
	global_options.add(command_options).add(config_simulation).add(config_field).add(config_plasma).add(config_impurities).add(config_particle);

	// Store variable map
	store(po::command_line_parser(ac,av).options(global_options).run(),
		vm);
	if(vm.count("help")){ // help messege
		po::options_description visible("Options description");
		if(vm.count("verbose")) visible.add(global_options);
		else visible.add(command_options);
		clog<< visible <<'\n';
		throw 0;
	}
	notify(vm);

	// Store input file
	std::ifstream ifs(indir.c_str()); // input file
	std::ifstream ofs(outdir.c_str()); // output directory
	if(!ifs) throw string("Cannot open file: ")+indir;
	else if(!ofs) throw string("Cannot open directory: ")+outdir;
	else{
	    store(parse_config_file(ifs,global_options),vm);
	    notify(vm);
	}

	if(vm.count("simulation_type")){ // Assign simulation type
		const string user_input= vm["simulation_type"].as<string>();
		if(user_input=="traj") swt.simu= SIM_TRAJ;
		else if(user_input=="ense") swt.simu= SIM_ENSE;
		else if(user_input=="poin") swt.simu= SIM_POIN;
		else throw string("Invalid simulation type: ")+user_input;
	}

	try{ // Store slices file
		if(vm.count("slices")){
			swt.flux= true;
			total_slices= // number of slice files
				(vm["slices"].as<vector<string>>()).size();
			std::clog <<"Reading slice files...\n";
			for(int i=0;i<total_slices;i++){
				const std::string slice_file=
					vm["slices"].as<vector<string>>()[i];
				std::ifstream ifsurf(slice_file.c_str());
				if(!ifsurf) throw "Cannot open file: "+slice_file;
			}
			read_slice();
		}
		else if(swt.simu==SIM_ENSE){
			throw string("Particle flux won't be calculated!");
		}
	}
	catch(const std::string& ex){
		clog <<"[WARNING] "+ex+"\n";
	}

	if(vm.count("plasma.temperature")){ // store plasma temperature
		const string temp_input=
			vm["plasma.temperature"].as<coordinate>().coord_first;
		if(temp_input=="flat") plasma.temp_prof= TEM_FLAT;
		else if(temp_input=="gaus") plasma.temp_prof= TEM_GAUS;
		else throw string("Invalid temperature profile: ") +temp_input;
		plasma.i_temp=
			stoi(vm["plasma.temperature"].as<coordinate>().coord_second);
	}
	if(vm.count("plasma.density")){ // store plasma density
		const string dens_input=
			vm["plasma.density"].as<coordinate>().coord_first;
		if(dens_input=="flat") plasma.dens_prof= DEN_FLAT;
		else if(dens_input=="para") plasma.dens_prof= DEN_PARA;
		else throw string("Invalid density profile: ") +dens_input;
		plasma.d_dens=
			stod(vm["plasma.density"].as<coordinate>().coord_second);
	}

	if(vm.count("density")) // store plasma density
		plasma.d_dens= vm["density"].as<double>();

	if(vm.count("impurities.type")){ // store impurities
		total_species=
			(vm["impurities.type"].as<vector<string>>()).size();
		for(int i=0;i<total_species;i++){
			const string species_type= // store species type
				vm["impurities.type"].as<vector<string>>()[i];
			if(species_type=="ions") species[i].sign= CHG_IONS;
			else if(species_type=="electrons") species[i].sign= CHG_ELEC;
			species[i].i_mass= // store species mass
				vm["impurities.mass"].as<vector<int>>()[i];
			species[i].i_chrg= // store species charge
				vm["impurities.charge"].as<vector<int>>()[i];

			const string temp_input= // temperature profile
				vm["impurities.temperature"].as<vector<coordinate>>()[i]
				.coord_first;
			if(temp_input=="flat") species[i].temp_prof= TEM_FLAT;
			else if(temp_input=="gaus") species[i].temp_prof= TEM_GAUS;
			else throw string("Invalid temperature profile: ") +temp_input;
			species[i].i_temp=
				stoi(vm["impurities.temperature"].as<vector<coordinate>>()[i]
					.coord_second);

			const string dens_input= // density profile
				vm["impurities.density"].as<vector<coordinate>>()[i]
				.coord_first;
			if(dens_input=="flat") species[i].dens_prof= DEN_FLAT;
			else if(dens_input=="para") species[i].dens_prof= DEN_PARA;
			else throw string("Invalid density profile: ") +dens_input;
			species[i].d_dens=
				stod(vm["impurities.density"].as<vector<coordinate>>()[i]
					.coord_second);
		}
	}

	// Assign values
	if(vm.count("collisions")) swt.coll= true;
	if(vm.count("perturbation")) swt.pert= true;
	if(vm.count("screen-log")) swt.scrn= true;
	if(vm.count("verbose")) swt.verb= true;
	if(vm.count("device")){ // device type
		const string user_input=
			vm["device"].as<string>();
		if(user_input=="tkmk") swt.conf= DEV_TKMK;
		else if(user_input=="stel") swt.conf= DEV_STEL;
		else throw string("device type: ")+user_input;
	}
	if(vm.count("simulation.poincare_bounds")){ // poincare bounds
		sim.d_rlow=
			stod(vm["simulation.poincare_bounds"]
			.as<coordinate>().coord_first);
		sim.d_rupp=
			stod(vm["simulation.poincare_bounds"]
			.as<coordinate>().coord_second);
	}
	if(vm.count("field.electrostatic_potential")){ // electrostatic potential
		const string user_input=
			vm["field.electrostatic_potential"].as<coordinate>().coord_first;
		if(user_input=="cons") field.potp= POT_CONS;
		else if(user_input=="pres") field.potp= POT_PRES;
		else throw string("Invalid potential profile: ")+user_input;
		field.d_epot=
			stod(vm["field.electrostatic_potential"]
			.as<coordinate>().coord_second);
	}
	if(vm.count("particle.type_test")){ // test particle type
		const string user_input=
			vm["particle.type_test"].as<string>();
		if(user_input=="ions") test.sign= CHG_IONS;
		else if(user_input=="electrons") test.sign= CHG_ELEC;
		else throw string("Invalid particle type: ")+user_input;
	}
	if(vm.count("particle.pitch_angle")){ // pitch distribution
		const string user_input=
			vm["particle.pitch_angle"].as<coordinate>().coord_first;
		if(user_input=="par") test.ptch_dist= PCH_PAR;
		else if(user_input=="ant") test.ptch_dist= PCH_ANT;
		else if(user_input=="iso") test.ptch_dist= PCH_ISO;
		else if(user_input=="fix"){
			test.ptch_dist=
				PCH_FIX;
			test.i_ptch=
				stoi(vm["particle.pitch_angle"].as<coordinate>().coord_second);
		}
		else throw string("Invalid pitch distribution: ")+user_input;
	}
	if(vm.count("particle.energy")){ // test particle energy
		const string user_input=
			vm["particle.energy"].as<coordinate>().coord_first;
		if(user_input=="mono") test.ener_dist= ENG_MONO;
		else if(user_input=="maxw") test.ener_dist= ENG_MAXW;
		else throw string("Invalid energy distribution: ")+user_input;
		test.d_ener=
			stod(vm["particle.energy"].as<coordinate>().coord_second);
	}

}


template <typename T> // bool condition exception
int ex_bool(const T& val, const std::string& messege) {
	if(val) throw std::invalid_argument(messege);
}

template <typename T, typename U> // upper bound exception
int ex_upper_bound(const T& val, const U& max, const std::string& messege){
	if(val>max){
		throw std::invalid_argument(
			"Invalid "+messege+": "+std::to_string(val)
			+". Max is: "+std::to_string(max));
	}
}

template <typename T, typename U> // lower bound exception
int ex_lower_bound(const T& val, const U& min, const std::string& messege){
	if(val<min){
		throw std::invalid_argument(
			"Invalid "+messege+": "+std::to_string(val)
			+". Min is: "+std::to_string(min));
	}
}

void check_bounds(){ // Bound input values
	ex_upper_bound(sim.i_part,MAX_PART,"number_of_particles");
	ex_upper_bound(sim.i_tran,MAX_TRAN,"toroidal_transits");
	ex_upper_bound(sim.i_sdel,MAX_SDEL,"step_delta");
	ex_upper_bound(sim.i_traj,MAX_TRAJPNTS,"trajectory_points");
	ex_upper_bound(sim.i_mome,MAX_MOMEPNTS,"moment_points");
	ex_upper_bound(sim.i_dist,MAX_DISTPNTS,"toroidal_transits");
	ex_upper_bound(sim.i_part,MAX_PART,"number_of_particles");
	ex_upper_bound(sim.i_tran,MAX_TRAN,"toroidal_transits");
	ex_upper_bound(sim.d_rupp,1.0,"poincare_upper_bound");
	ex_lower_bound(sim.d_rlow,0.0,"poincare_lower_bound");
	ex_bool(sim.d_rupp<sim.d_rlow,
		"poincare_upper_bound should be larger than poincare_lower_bound");
	ex_lower_bound(field.i_baxi,0.0,"field_magnitude");
	ex_lower_bound(field.i_qaxi,0,"q_on_axis");
	ex_lower_bound(field.i_qwll,0,"q_on_wall");
	ex_bool(field.i_qwll<field.i_qaxi,
		"q_on_wall should be larger than q_on_axis");
	ex_lower_bound(field.i_rmin,0,"minor_radius");
	ex_lower_bound(field.i_rmaj,0,"major_radius");
	ex_bool(field.i_rmaj<field.i_rmin,
		"major_radius should be larger than minor_radius");
	ex_lower_bound(field.d_ampx,0.0,"perturbation_amplitude");
	ex_lower_bound(field.i_freq,0,"perturbation_phase");
	ex_lower_bound(field.i_mmod,0,"m_mode");
	ex_lower_bound(field.i_nmod,0,"n_mode");
	// species
	ex_upper_bound(total_species,MAX_SPECIES,"total_species");
	for(int i=0;i<total_species;i++)
	{
		ex_lower_bound(species[i].i_mass,0,"mass");
		ex_lower_bound(species[i].i_chrg,0,"charge");
		ex_lower_bound(species[i].i_temp,0,"temperature");
		ex_lower_bound(species[i].d_dens,0.0,"density");
	}

	ex_lower_bound(test.i_mass,0,"test_mass");
	ex_lower_bound(test.i_chrg,0,"test_charge");
	ex_upper_bound(test.i_ptch,360,"pitch_angle");
	ex_lower_bound(test.i_ptch,0,"pitch_angle");
	ex_lower_bound(test.d_ener,0.0,"energy");
	ex_upper_bound(test.d_radi,1.0,"radial_position");
	ex_lower_bound(test.d_radi,0.0,"radial_position");
	ex_upper_bound(test.i_thet,180,"theta");
	ex_lower_bound(test.i_thet,-180,"theta");
	ex_upper_bound(test.i_zeta,360,"zeta");
	ex_lower_bound(test.i_zeta,0,"zeta");
}


std::fstream dlog; // output variable
template <typename T>
void LOG_VAL( // logger function
	const std::string &sec,
	const std::string& str,
	const T& val){
	dlog <<std::left <<std::setw(24) <<str <<" = ";
	dlog <<std::left <<std::setw(12)
		<<std::setprecision(2) <<std::scientific <<val;
	const std::string item= sec+"."+str;
	if(vm.count(item)) dlog <<std::right <<std::setw(10) <<"# user\n";
	else dlog <<std::setw(12) <<std::right <<"# default\n";
	if(swt.scrn){
		std::clog <<std::left <<std::setw(24) <<str <<" = ";
		std::clog <<std::left <<std::setw(12)
			<<std::setprecision(2) <<std::scientific <<val;
		if(vm.count(item)) std::clog <<std::right <<std::setw(10) <<"# user\n";
		else std::clog <<std::setw(12) <<std::right <<"# default\n";
	}
}
void LOG_VAL(
	const std::string& str,
	const std::string& val){
	dlog <<std::left <<std::setw(24) <<str <<" = ";
	dlog <<std::left <<std::setw(8) <<std::setprecision(2)
		 <<std::scientific <<val <<'\n';
	if(swt.scrn){
		std::clog <<std::left <<std::setw(24) <<str <<" = ";
		std::clog <<std::left <<std::setw(8) <<std::setprecision(2)
			 <<std::scientific <<val <<'\n';
	}
}

void LOG_STR(const std::string& str){
	if(swt.scrn) std::clog <<std::left <<str <<'\n';
	dlog <<std::left <<str <<'\n';
}

void logger(){ // record log
	namespace po = boost::program_options;
	std::ofstream(outdir+"/data.log");
	dlog.open(outdir+"/data.log",std::ios::app);
	LOG_STR("# Log");

	// for (const auto& it : vm) {
	// 	std::string first= it.first.c_str();
	// 	auto& value = it.second.value();
	// 	if(auto v= boost::any_cast<int>(&value))
	// 		LOG_VAL(first,*v);
	// 	else if(auto v= boost::any_cast<double>(&value))
	// 		LOG_VAL(first,*v);
	// 	else if(auto v= boost::any_cast<coordinate>(&value))
	// 		{LOG_VAL(first,(*v).coord_first);}
	// 	else if(auto v= boost::any_cast<std::string>(&value)){
	// 		if(*v=="") LOG_STR("[USER]    "+first);
	// 		else LOG_VAL(first,*v);
	// 		}
	// 	else
	// 		throw std::string("Invalid parsing type");
	// }
	LOG_STR("\n[switches]");
	switch(swt.simu){
		case SIM_TRAJ: LOG_VAL("simulation_type","traj"); break;
		case SIM_ENSE: LOG_VAL("simulation_type","ense"); break;
		case SIM_POIN: LOG_VAL("simulation_type","poin"); break;
	}
	if(swt.coll) LOG_VAL("collisions","on");
	else LOG_VAL("collisions","off");
	if(swt.pert) LOG_VAL("perturbation","on");
	else LOG_VAL("perturbation","off");
	switch (swt.conf) {
		case DEV_TKMK: LOG_VAL("device","tkmk"); break;
		case DEV_STEL: LOG_VAL("device","stel"); break;
	}

	LOG_STR("\n[simulation]");
	LOG_VAL("simulation","number_of_particles",sim.i_part);
	LOG_VAL("simulation","toroidal_transits",sim.i_tran);
	LOG_VAL("simulation","step_delta",sim.i_sdel);
	LOG_VAL("simulation","trajectory_points",sim.i_traj);
	LOG_VAL("simulation","moment_points",sim.i_mome);
	LOG_VAL("simulation","distribution_points",sim.i_dist);
	LOG_VAL("simulation","poincare_bounds",
		vm["simulation.poincare_bounds"].as<coordinate>().coord_first +" "
		+vm["simulation.poincare_bounds"].as<coordinate>().coord_second);
	LOG_VAL("simulation","poincare_plane",sim.i_surf);

	LOG_STR("\n[field]");
	LOG_VAL("field","field_magnitude",field.i_baxi);
	LOG_VAL("field","q_on_axis",field.i_qaxi);
	LOG_VAL("field","q_on_wall",field.i_qwll);
	LOG_VAL("field","major_radius",field.i_rmaj);
	LOG_VAL("field","minor_radius",field.i_rmin);
	LOG_VAL("field","perturbation_amplitude",field.d_ampx);
	LOG_VAL("field","perturbation_frequency",field.i_freq);
	LOG_VAL("field","perturbation_phase",field.i_phas);
	LOG_VAL("field","m_mode",field.i_mmod);
	LOG_VAL("field","n_mode",field.i_nmod);
	LOG_VAL("field","electrostatic_potential",
		vm["field.electrostatic_potential"].as<coordinate>().coord_first +" "
		+vm["field.electrostatic_potential"].as<coordinate>().coord_second);

	LOG_STR("\n[impurities]");
	for(int i=0;i<total_species;i++)
	{
		LOG_VAL("impurities","# type",
			vm["impurities.type"].as<std::vector<std::string>>()[i]);
		LOG_VAL("impurities","# mass",species[i].i_mass);
		LOG_VAL("impurities","# charge",species[i].i_chrg);
		LOG_VAL("impurities","# temperature",
			vm["impurities.temperature"].as<std::vector<coordinate>>()[i]
			.coord_first+" "
			+vm["impurities.temperature"].as<std::vector<coordinate>>()[i]
			.coord_second);
		std::ostringstream streamObj;
		streamObj << species[i].d_dens;
		std::string strObj = streamObj.str();
		LOG_VAL("impurities","# density",
			vm["impurities.density"].as<std::vector<coordinate>>()[i]
			.coord_first+" "
			+strObj);
	}

	LOG_STR("\n[particle]");
	LOG_VAL("particle","type_test",
		vm["particle.type_test"].as<std::string>());
	LOG_VAL("particle","mass_test",test.i_mass);
	LOG_VAL("particle","charge_test",test.i_chrg);
	LOG_VAL("particle","radial_position",test.d_radi);
	LOG_VAL("particle","theta",test.i_thet);
	LOG_VAL("particle","zeta",test.i_zeta);
	LOG_VAL("particle","pitch_angle",
		vm["particle.pitch_angle"].as<coordinate>().coord_first +" "
		+vm["particle.pitch_angle"].as<coordinate>().coord_second);
	LOG_VAL("particle","energy",
		vm["particle.energy"].as<coordinate>().coord_first +" "
		+vm["particle.energy"].as<coordinate>().coord_second);
	dlog.close();
}


void info(){ //print info
	const char* myinfo =
R"(gcafpp code by Rodrigo Saavedra
Email: saavestro@gmail.com)";
	std::clog <<myinfo <<"\n";
}
