//
//  pde_solver.hpp
//  PDE_solver
//
//  Created by Florian on 08/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#ifndef pde_solver_hpp
#define pde_solver_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "Option.hpp"
#include "Payoff.h"
#include "pde.hpp"



namespace Solve
{

    class BS_Solver
    {
	public:
		//BS_Solver();
		BS_Solver(BS_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity);
		void calculate_parameters();
		void set_initial_conditions();
		virtual std::vector<double> boundary_increment(const double& t);

        virtual std::vector<double> forward_coefficient(const double& temp, const size_t& i = 0);
        virtual std::vector<double> present_coefficient(const double& temp, const size_t& i = 0);
        virtual std::vector<double> backward_coefficient(const double& temp, const size_t& i = 0);
		dauphine::matrix transition_matrix(const double& temp, const size_t& i = 0);
		virtual void Crout_Algo_Resolution();
		std::vector<double> LU_compute( dauphine::matrix& L, dauphine::matrix& U, const std::vector<double>& b);
		std::vector<double> get_option_payoff();
		std::vector<double> get_S_grid();
        std::vector<double> get_price_curve();
		double get_price(const double& S);
		double compute_delta(const double& S);
		std::vector<double> compute_delta();
		double compute_gamma(const double& S);
		std::vector<double> compute_gamma();
		double compute_theta(const double& S);
		std::vector<double> compute_theta();
		std::vector<double> compute_vega();
		double compute_vega(const double& S);

	private: 
		BS_PDE* pde;
		double theta;
		double S0;
		double x_max;
		double x_min;
		double dx;
		std::size_t space_dim;
		std::vector<double> x_values;
		std::vector<double> S_values;
		double dt;
		double maturity;
		std::size_t time_dim;
		std::string l;
		std::string r;
		std::vector<double> option_payoff;
		std::vector<double> old_result;
		std::vector<double> new_result;
		bool resolved;
		std::vector<double> vega;
    };

	class Exo_Solver : public BS_Solver
	{
	public:
		Exo_Solver(Exo_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity);
		std::vector<double> forward_coefficient(const double& temp, const size_t& i);
		std::vector<double> present_coefficient(const double& temp, const size_t& i);
		std::vector<double> backward_coefficient(const double& temp, const size_t& i);
		std::vector<double> boundary_increment(const size_t& i);
		void Crout_Algo_Resolution();
	private:
		Exo_PDE* pde;
		double theta;
		double S0;
		double x_max;
		double x_min;
		double dx;
		std::size_t space_dim;
		std::vector<double> x_values;
		std::vector<double> S_values;
		double dt;
		double maturity;
		std::size_t time_dim;
		std::string l;
		std::string r;
		std::vector<double> option_payoff;
		std::vector<double> old_result;
		std::vector<double> new_result;
		bool resolved;
		//std::vector<double> vega;
	};


}
#endif /* pde_solver_hpp */
