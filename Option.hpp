#pragma once


#ifndef Option_hpp
#define Option_hpp

#include "Payoff.h"
#include "matrix.hpp"
#include <vector>

class VanillaOption {
	public:
		PayOff* pay_off;

		const double K = 0.0;
		const double r = 0.0;
		const double T = 0.0;
		const double sigma = 0.0;



		VanillaOption(const double& _K, const double& _r, const double& _T,
			const double& _sigma, PayOff* _pay_off);
		VanillaOption* Option_vega( const double& dv);
};

// We wanted to use inheritance here too but had some issues with the need to define a default constructor, non access to attributes.

class ExoticOption {
	private:
		std::vector<double> R;
		std::vector<double> Sigma;

	public:
		PayOff* pay_off;

		const double K = 0.0;
		const double r = 0.0;
		const double T = 0.0;

		ExoticOption(const double& _K, const std::vector<double>& _r, const double& _T,
			const std::vector<double>& _sigma, PayOff* _pay_off);
		std::vector<double> get_yield_curve();
		std::vector<double> get_vol_TS();
};
#endif
