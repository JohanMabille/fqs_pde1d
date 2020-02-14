#ifndef __VANILLA_OPTION_CPP
#define __VANILLA_OPTION_CPP

#include "Option.hpp"



VanillaOption::VanillaOption(const double& _K, const double& _r, const double& _T,
	const double& _sigma, PayOff* _pay_off) :
	K(_K), r(_r), T(_T), sigma(_sigma), pay_off(_pay_off) {}


VanillaOption* VanillaOption::Option_vega(const double& dv) 
{
	VanillaOption* vegaoption = new VanillaOption(K, r, T, sigma + dv, pay_off);
	return vegaoption;
}

ExoticOption::ExoticOption(const double& _K, const std::vector<double>& _r, const double& _T,
	const std::vector<double>& _sigma, PayOff* _pay_off) :
	VanillaOption(_K, 0.0, _T, 0.0, _pay_off), Sigma(_sigma), R(_r), pay_off(_pay_off), K(_K),
	T(_T){}

std::vector<double> ExoticOption::get_yield_curve(){
	return R;
}

std::vector<double> ExoticOption::get_vol_TS(){
	return Sigma;
}



#endif
