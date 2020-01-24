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
#include <fstream>
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

    std::ofstream outFile;
    outFile.open("PDE_Solver_Outputs.csv");
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return -1;
    }

    std::vector<double> S_grid = PDE_solve->get_S_grid();
    std::vector<double> price_curve = PDE_solve->get_price_curve();
    std::vector<double> delta_curve = PDE_solve->compute_delta();
    std::vector<double> gamma_curve = PDE_solve->compute_gamma();
    std::vector<double> theta_curve = PDE_solve->compute_theta();
    std::vector<double> vega_curve = PDE_solve->compute_vega();
    std::vector<double> payoff_curve = PDE_solve->get_option_payoff();

   
    outFile << "Spot,Option Price,Delta,Gamma,Theta,Vega" << std::endl;
    for (size_t i = 0; i < space_dim - 1; ++i) 
    {
        outFile << S_grid[i] << ","
            << price_curve[i] << "," 
            << delta_curve[i] << ","
            <<gamma_curve[i] << ","
            <<theta_curve[i] << ","
            <<vega_curve[i]<< ","
            <<payoff_curve[i]
            << std::endl;
    }

    outFile << delta_curve[80] << std::endl;
    

    

    outFile.close();

    return 0;
    
}
