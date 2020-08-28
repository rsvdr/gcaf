//----------------------------------------------------------------//
// gcafpp code by R. Saavedra. Email: saavestro@gmail.com         //
// For more information visit https://github.com/Saavestro/gcafpp //
//----------------------------------------------------------------//
#include "shared.h"
#include <random>
#include <algorithm>
#include <iterator>

const double urand();

const double normal_dist(
	const double& a,
	const double& b){ //normal distribution
	const double rone= urand();
	const double rtwo= urand();
	const double x= sqrt(-2*log(rone))*cos(2*M_PI*rtwo);
	return a +b*x;
}

const double urand(){ // random number from 0.0 to 1.0
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> frand(0.0, 1.0);
	return frand(mt);
}

const double signnum(const double& x){ //sign function
  if (x>0.0) return 1.0;
  if (x<0.0) return -1.0;
  return x;
}

const double xproj(
	const double& px,
	const double& tx,
	const double& zx){ //x projection
	const double pdum= 2.0*code.d_leng*code.d_leng*px;
	const double rdum= sqrt(pdum);
	return (1.0*field.i_rmaj +rdum*cos(tx))*cos(zx);
}

const double yproj(
	const double& px,
	const double& tx,
	const double& zx){ //y projection
	const double pdum= 2.0*code.d_leng*code.d_leng*px;
	const double rdum= sqrt(pdum);
	return (1.0*field.i_rmaj +rdum*cos(tx))*sin(zx);
}

const double zproj(
	const double& px,
	const double& tx,
	const double& zx){ //z projection
	const double pdum= 2.0*code.d_leng*code.d_leng*px;
	const double rdum= sqrt(pdum);
	return rdum*sin(tx);
}

std::vector<std::string> strip(
	std::string lngarg){ //string elements into array
	std::istringstream iss(lngarg);
	std::vector<std::string> tokens;
	std::copy(
		std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter(tokens)
	);
	if(tokens.size()==0) tokens.push_back("EMPTY");
	return tokens;
}
