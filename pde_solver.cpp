//
//  pde_solver.cpp
//  PDE_solver
//
//  Created by Florian on 08/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "pde_solver.hpp"

namespace Solve
{

    BS_Solver::BS_Solver(BS_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity)
    :pde(_pde), theta(_theta), space_dim(_space_dim), time_dim(_time_dim), S0(_S0), maturity(_maturity) {
		BS_Solver::calculate_parameters();
	}
 
	void BS_Solver::calculate_parameters() {
		double stdv = pde->standard_dev();

                // Design: this could be in a dedicated class (mesh, or grid)
		x_max = log(S0) + 5 * stdv;
		x_min = log(S0) - 5 * stdv;
		dx = 10 * stdv / static_cast<double>(space_dim);
		dt =  maturity / static_cast<double>(time_dim);
		r = pde->get_right_boundary_type();
		l = pde->get_left_boundary_type();
		resolved = false;

		x_values.resize(space_dim - 1, 0.0);
		S_values.resize(space_dim - 1, 0.0);
		double val = x_min;
		//filling x_values ;
		for (auto x = x_values.begin(); x != x_values.end(); ++x)
		{
			*x = val + dx;
			val = *x;
		}
		//filling S_values
		std::transform(x_values.begin(), x_values.end(), S_values.begin(),
			[](double d) { return exp(d); });

	}
    
    std::vector<double> BS_Solver::forward_coefficient(const double& temp, const size_t& i)
    {
        // Implementation: double dx_2 = dx * dx is way faster
        double dx_2=pow(dx, 2.0);
        std::vector<double>a_coef(space_dim-2);
        
        for(auto it = a_coef.begin(); it != a_coef.end(); ++it)
        {
			*it = temp *(-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);             
        }
		return a_coef;
	}
    
    std::vector<double> BS_Solver::present_coefficient(const double& temp, const size_t& i)
    {
        double dx_2=pow(dx,2.0);
        std::vector<double>b_coef(space_dim-1);
		//auto ptr = b_coef.begin();
		//auto ftr = b_coef.back();

                // Implementation: For efficiency, prefer comparing enum values than string
                // Design: consider providing an hierarchy of classes for boundary conditions.
                // The user may want to combine them (Dirichlet on left wing, Neuman on right wing for instances)
                // Also your handling of Boundary conditions is way too complicated. Compute your tridiag system
                // M X = B for i in [1, S - 2], then directly fill m(0, 0), m(0, 1), m(S-1, S-2), m(S-1, S-1) and
                // X[0] and X[1] with boundary conditions. For Dirichlet this is simply:
                // m(0, 0) = 1; m(0, 1) = 0; X[0] = previous_x[0];
                // m(S-1, S-2) = 0; m(S-1, S-1) = 1; X[S-1] = previous_x[S-1]
		if (l.compare("N") == 0) 
		{
			//present coeff + backward coeff because of neumann condition
			b_coef[0] = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2) + temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		else
		{
			b_coef[0] = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
		}

		for (std::size_t i = 1; i < space_dim - 2; ++i)
		{
			b_coef[i] = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
		}

		if (r.compare("N") == 0)
		{
			//present coeff + forward coeff because of neumann condition
			b_coef[space_dim - 2] = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2) + temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		else
		{
			b_coef[space_dim - 2] = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
		}

		return b_coef;
	}
    
    std::vector<double> BS_Solver::backward_coefficient(const double& temp, const size_t& i)
    {
		double dx_2 = pow(dx, 2.0);
		std::vector<double> c_coef(space_dim - 2);

		for (auto it = c_coef.begin(); it != c_coef.end(); ++it)
		{
			*it = temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		return c_coef;
    }
    
	void BS_Solver::set_initial_conditions() 
	{
		
		option_payoff.resize(space_dim - 1, 0.0);
		//filling  and initial condition in result vector
		
		option_payoff = pde->init_cond(S_values);
	}

	std::vector<double> BS_Solver::boundary_increment(const double& t)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double> boundary_effect;
		boundary_effect.resize(space_dim -1, 0.0);
		//auto it1 = boundary_effect.begin(); //check this with prof
		//auto it2 = boundary_effect.back();
		double left_b = pde->boundary_left();
		double right_b = pde->boundary_right(t, exp(x_max));

		if (l.compare("D") == 0)
		{
			boundary_effect[0] = - left_b * dt * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the backward coefficient.
		}
		else
		{
			boundary_effect[0] = left_b * dx * dt * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the backward coefficient.
		}

		if (r.compare("D") == 0)
		{
			boundary_effect[space_dim - 2] = -right_b * dt * (-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the forward coefficient.
		}
		else
		{
			boundary_effect[space_dim - 2] = right_b * dx * dt * (-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the forward coefficient.
		}
		return boundary_effect;
	}
    
	dauphine::matrix BS_Solver::transition_matrix(const double& temp, const size_t& i)
	{
                // Implementation: why rvalue references here?
		std::vector<double>&& a = forward_coefficient(temp, i);
		std::vector<double>&& b = present_coefficient(temp, i);
		std::vector<double>&& c = backward_coefficient(temp, i);
		dauphine::matrix m(space_dim - 1, space_dim - 1);

		m(0, 0) = b[0];
		m(0, 1) = a[0];

		for (std::size_t i = 1; i != space_dim - 2 ; ++i)
		{
			m(i, i - 1) = c[i-1];
			m(i, i) = b[i];
			m(i, i + 1) = a[i];
		}

		m(space_dim - 2, space_dim - 3) = c[space_dim - 3];
		m(space_dim - 2, space_dim - 2) = b[space_dim - 2];

		return m ;
	}

	void BS_Solver::Crout_Algo_Resolution()
	{
		double temp_lhs = theta * dt;
		double temp_rhs = - (1- theta) * dt;

		//setting initial conditions
		set_initial_conditions();


		//creation of the left and right transition matrices
                // Implementation: Why rvalue references here?
                // Implementation: consider allocating your matrices
                // and vectors once and for good in this main function,
                // and pass them by references to other routines so they can
                // fill them. This would avoid a lot of allocations and copies
		dauphine::matrix&& M_lhs = transition_matrix(temp_lhs);
		dauphine::matrix&& M_rhs = transition_matrix(temp_rhs);

		// L - U decomposition of the left matrix for the Crout Algorithm 
		dauphine::matrix L(space_dim -1, space_dim - 1);
		dauphine::matrix U(space_dim - 1, space_dim - 1);
		L(0, 0) = M_lhs(0, 0);
		U(0, 0) = 1.0;
		U(0, 1) = M_lhs(0, 1) / L(0, 0);

		for (std::size_t i = 1; i != space_dim - 2; ++i)
		{
			L(i, i - 1) = M_lhs(i, i - 1);
			L(i, i) = M_lhs(i, i) - L(i, i - 1) * U(i - 1, i);
			U(i, i + 1) = M_lhs(i, i + 1) / L(i, i);
			U(i, i) = 1.0;
		}

		U(space_dim - 2, space_dim - 2) = 1.0;
		L(space_dim - 2, space_dim - 3) = M_lhs(space_dim - 2, space_dim - 3);
		L(space_dim - 2, space_dim - 2) = M_lhs(space_dim - 2, space_dim - 2) - L(space_dim - 2, space_dim - 3) * U(space_dim - 3, space_dim - 2);
		
		std::vector<double> v(space_dim - 1);
		std::vector<double> tmp(space_dim - 1);
		double t;
		new_result = option_payoff;

		for (std::size_t i = 1; i != time_dim; ++i)
		{
			old_result = new_result;
			t = maturity - i * dt; // think about && rvalue lvalue;
			tmp = boundary_increment(t);
			v = M_rhs.produit_mat_vect(old_result);
			std::transform(tmp.begin(), tmp.end(), v.begin(), v.begin(), std::plus<double>());

			new_result = LU_compute(L, U, v);
		}
		resolved = true;
	}

	std::vector<double> BS_Solver::LU_compute(dauphine::matrix& L, dauphine::matrix& U, const std::vector<double>& b)
	{
		std::vector<double> x(b.size(),0.0);
		std::vector<double> z(b.size(),0.0);
		
		z[0] = b[0] / L(0, 0);

		for (std::size_t i = 1; i != b.size(); ++i)
		{
			z[i] = (b[i] - L(i, i - 1) * z[i - 1]) / L(i, i);
		}

		x[b.size() - 1] = z[b.size() - 1];
		for (auto i = b.size() - 2; i != 0; --i)
		{
			x[i] = z[i] - U(i, i + 1) * x[i + 1];
		}

		return x;
	}

	std::vector<double> BS_Solver::get_option_payoff()
	{
		return option_payoff;
	}

	std::vector<double> BS_Solver::get_S_grid()
	{
		return S_values;
	}

	std::vector<double> BS_Solver::get_price_curve()
	{
		//test if resolution of the PDE has been done
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		return new_result;
	}

	double BS_Solver::get_price(const double& S)
	{
		//add test if resolution of the PDE has been done
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		double x = log(S);
		int i = 0;
		while (x_values[i] < x)
		{
			++i;
		}
		return new_result[i];
	}

	double BS_Solver::compute_delta(const double& S)
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		int i = 0;
		while (S_values[i] < S)
		{
			++i;
		}
		if (i == 0) {
			return (new_result[1] - new_result[0]) / (S_values[1] - S_values[0]);
		}
		else if (i == space_dim - 2) {
			return (new_result[space_dim - 2] - new_result[space_dim - 3]) / (S_values[space_dim - 2] - S_values[space_dim - 3]);
		}
		else {
			return (new_result[i + 1] - new_result[i - 1]) / (S_values[i + 1] - S_values[i - 1]);
		}
	}

	std::vector<double> BS_Solver::compute_delta()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> delta;
		for (size_t i = 0; i < space_dim - 1;++i) 
		{
			delta.push_back(compute_delta(S_values[i]));
		}
		
		
		return delta;
	}

	double BS_Solver::compute_gamma(const double& S)
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		int i = 0;
		while (S_values[i] < S)
		{
			++i;
		}
		//Discretization of second derivative here is different because S-grid is non uniform.
		double ds_2;
		if (i == 0) {
			ds_2 = (S_values[2] - S_values[1]) * (S_values[1] - S_values[0]) * (S_values[2] - S_values[0]);
			return 2 / ds_2 * ((new_result[2] - new_result[1]) * (S_values[1] - S_values[0]) - (new_result[1] - new_result[0]) * (S_values[2] - S_values[1]));
		}
		else if (i == space_dim - 2) {
			ds_2 = (S_values[space_dim - 2] - S_values[space_dim - 3]) * (S_values[space_dim -3] - S_values[space_dim - 4]) * (S_values[space_dim - 2] - S_values[space_dim - 4]);
			return 2 / ds_2 * ((new_result[space_dim - 2] - new_result[space_dim - 3]) * (S_values[space_dim - 3] - S_values[space_dim - 4]) - (new_result[space_dim - 3] - new_result[space_dim - 4]) * (S_values[space_dim - 2] - S_values[space_dim - 3]));
		}
		else {
			ds_2 = (S_values[i + 1] - S_values[i]) * (S_values[i] - S_values[i - 1]) * (S_values[i + 1] - S_values[i - 1]);
			return 2 / ds_2 * ((new_result[i + 1] - new_result[i]) * (S_values[i] - S_values[i - 1]) - (new_result[i] - new_result[i - 1]) * (S_values[i + 1] - S_values[i]));
		}
	}

	std::vector<double> BS_Solver::compute_gamma()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> gamma;
		for (size_t i = 0; i < space_dim - 1; ++i)
		{
			gamma.push_back(compute_gamma(S_values[i]));
		}
		return gamma;
	}

	double BS_Solver::compute_theta(const double& S)
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		double x = log(S);
		int i = 0;
		while (x_values[i] < x)
		{
			++i;
		}
		return (old_result[i] - new_result[i]) * dt * 365; //theta per 1 day change
	}

	std::vector<double> BS_Solver::compute_theta()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> theta;
		for (size_t i = 0; i < space_dim - 1; ++i)
		{
			theta.push_back(compute_theta(S_values[i]));
		}
		return theta;
	}
	
	std::vector<double> BS_Solver::compute_vega()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		vega.resize(space_dim - 1);
                // Implementation: memory leak, vega_pde1 is never deleted
		BS_PDE* vega_pde1 = pde->vega_pde(0.01);
                // Implementation: memory leak, vega_pde2 is never deleted
		BS_PDE* vega_pde2 = pde->vega_pde(-0.01);
                // Implementation: memory leak
		BS_Solver* vega_solve1 = new Solve::BS_Solver(vega_pde1, theta, space_dim, time_dim, S0, maturity);
                // Implementation: memory leak
		BS_Solver* vega_solve2 = new Solve::BS_Solver(vega_pde2, theta, space_dim, time_dim, S0, maturity);
		vega_solve1->Crout_Algo_Resolution();
		vega_solve2->Crout_Algo_Resolution();
		std::vector<double> vega_price1 = vega_solve1->get_price_curve();
		std::vector<double> vega_price2 = vega_solve2->get_price_curve();
		for (size_t i = 0; i < space_dim - 1; ++i) {
			vega[i] = ((vega_price1[i] - vega_price2[i]) / 2);
		}
		return vega;
	}

	double BS_Solver::compute_vega(const double& S)
	{
		if (vega.size() == 0){
			vega = compute_vega();
		}
		
		double x = log(S);
		int i = 0;
		while (x_values[i] < x)
		{
			++i;
		}
		return vega[i];
	}
	


        
    

}



