//
//  main.cpp
//  Project C++
//
//  Created by Florian on 02/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
//#include "Rates.hpp"
#include "Payoff.h"
#include "pde.hpp"
#include "pde_solver.hpp"
#include "matrix.hpp"



int main(int argc, char* argv[])
{
	double S0 = 100.0;
	double K = 100.0;
	double sigma = 0.15;
	double maturity = 1;
	double r = 0.01;
	double theta = 0.5;
    
    
    std::size_t space_dim = 100;
	std::size_t time_dim = 50;
	std::string l_boundary_type = "D";
	std::string r_boundary_type = "D";
	PayOff* pay_off_call = new PayOffCall(K);
	VanillaOption* call_option = new VanillaOption(K, r, maturity, sigma, pay_off_call);
	BS_PDE* bs_pde = new BS_PDE(call_option, l_boundary_type, r_boundary_type);
	Solve::matrix_pde_case1* PDE_solve = new Solve::matrix_pde_case1(bs_pde, theta, space_dim, time_dim, S0, maturity);
	PDE_solve->Crout_Algo_Resolution();
	double price = PDE_solve->get_price(S0);

    std::vector<double> price_curve = PDE_solve->get_price_curve();
    std::vector<double> delta_curve = PDE_solve->compute_delta();
    std::vector<double> gamma_curve = PDE_solve->compute_gamma();
    std::vector<double> theta_curve = PDE_solve->compute_theta();
    std::vector<double> vega_curve = PDE_solve->compute_vega();

    std::cout << "the price vector of the option at t=0 is :" << std::endl;
    for (auto it = price_curve.begin(); it != price_curve.end(); it++){
        std::cout << *it << std::endl;
    }

    std::cout << "the delta vector of the option at t=0 is :" << std::endl;
    for (auto it = delta_curve.begin(); it != delta_curve.end(); it++) {
        std::cout << *it << std::endl;
    }

    std::cout << "the gamma vector of the option at t=0 is :" << std::endl;
    for (auto it = gamma_curve.begin(); it != gamma_curve.end(); it++) {
        std::cout << *it << std::endl;
    }

    std::cout << "the theta vector of the option at t=0 is :" << std::endl;
    for (auto it = theta_curve.begin(); it != theta_curve.end(); it++) {
        std::cout << *it << std::endl;
    }

    std::cout << "the vega vector of the option at t=0 is :" << std::endl;
    for (auto it = vega_curve.begin(); it != vega_curve.end(); it++) {
        std::cout << *it << std::endl;
    }

    
    return 0;
    
}
