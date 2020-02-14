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

    // Design: the solver should be agnostic to the kind of options it is pricing
    // Actually it should even be agnostic to the "financial" aspect of the problem,
    // and simply solve the generalized heat equation
    // a(x, t) d2f/dx2 + b(x, t) df/dx + c(x, t) f + d(x, t)
    // You abstracted these coefficients in your PDE hierarchy, which is a good idea,
    // rely on this abstraction: a single PDE_solver that accepts a basicPDE argument
    // (which can actually be any kind or inheriting class, but the solver does not
    // need to know).
    class BS_Solver
    {
	public:
		BS_Solver(BS_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity);
                // Design: where is the virtual destructor?
                // Design: entity semantic
		virtual void calculate_parameters();
		virtual void set_initial_conditions();
		virtual std::vector<double> boundary_increment(const double& t);

                // implementation: consider passing builtin types by value
        virtual std::vector<double> forward_coefficient(const double& temp, const size_t& i = 0);
        virtual std::vector<double> present_coefficient(const double& temp, const size_t& i = 0);
        virtual std::vector<double> backward_coefficient(const double& temp, const size_t& i = 0);
		virtual dauphine::matrix transition_matrix(const double& temp, const size_t& i = 0);
		virtual void Crout_Algo_Resolution();
		virtual std::vector<double> LU_compute( dauphine::matrix& L, dauphine::matrix& U, const std::vector<double>& b);
		virtual std::vector<double> get_option_payoff();
		virtual std::vector<double> get_S_grid();
        virtual std::vector<double> get_price_curve();
		virtual double get_price(const double& S);
		virtual double compute_delta(const double& S);
		virtual std::vector<double> compute_delta();
		virtual double compute_gamma(const double& S);
		virtual std::vector<double> compute_gamma();
		virtual double compute_theta(const double& S);
		virtual std::vector<double> compute_theta();
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
                // Design: you're not on twitter, you can choose long and explicit names
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
		void calculate_parameters();
		void set_initial_conditions();
		std::vector<double> forward_coefficient(const double& temp, const size_t& i);
		std::vector<double> present_coefficient(const double& temp, const size_t& i);
		std::vector<double> backward_coefficient(const double& temp, const size_t& i);
		std::vector<double> boundary_increment(const size_t& i);
		void Crout_Algo_Resolution();

		virtual std::vector<double> get_price_curve();
		virtual double get_price(const double& S);
		virtual double compute_delta(const double& S);
		virtual std::vector<double> compute_delta();
		virtual double compute_gamma(const double& S);
		virtual std::vector<double> compute_gamma();
		virtual double compute_theta(const double& S);
		virtual std::vector<double> compute_theta();
		virtual std::vector<double> get_option_payoff();
	private:
                // Design: why duplicating data members if you use inheritance?
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
		std::vector<double> vega;
	};


}
#endif /* pde_solver_hpp */
