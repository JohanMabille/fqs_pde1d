#ifndef __VANILLA_OPTION_CPP
#define __VANILLA_OPTION_CPP

#include "option.hpp"

VanillaOption::VanillaOption() {}

VanillaOption::VanillaOption(double _K, double _r, double _T,
	double _sigma, PayOff* _pay_off) :
	K(_K), r(_r), T(_T), sigma(_sigma), pay_off(_pay_off) {}


VanillaOption* VanillaOption::Option_vega() 
{
	VanillaOption* vegaoption = new VanillaOption(K, r, T, sigma + 0.1, pay_off);
	return vegaoption;
}
#endif