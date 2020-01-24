#pragma once


#ifndef Option_hpp
#define Option_hpp

#include "Payoff.h"

class VanillaOption {
	public:
		PayOff* pay_off;

		const double K = 0.0;
		const double r = 0.0;
		const double T = 0.0;
		const double sigma = 0.0;



		VanillaOption(double _K, double _r, double _T,
			double _sigma, PayOff* _pay_off);
		VanillaOption* Option_vega( const double& dv);
};

#endif
