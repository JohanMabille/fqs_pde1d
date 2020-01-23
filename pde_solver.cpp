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
    matrix_pde_case1::matrix_pde_case1(BS_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity)
    :pde(_pde), theta(_theta), space_dim(_space_dim), time_dim(_time_dim), S0(_S0), maturity(_maturity) {
		calculate_parameters();
	}
 
	void matrix_pde_case1::calculate_parameters() {
		double stdv = pde->standard_dev();

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
			[](double d) { return std::exp(d); });

	}
    
    std::vector<double> matrix_pde_case1::forward_coefficient(const double& temp)
    {
        double dx_2=pow(dx, 2.0);
        std::vector<double>a_coef(space_dim-2);
        
        for(auto it = a_coef.begin(); it != a_coef.end(); ++it)
        {
			*it = temp *(-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);             
        }
		return a_coef;
	}
    
    std::vector<double> matrix_pde_case1::present_coefficient(const double& temp)
    {
        double dx_2=pow(dx,2.0);
        std::vector<double>b_coef(space_dim-1);
		//auto ptr = b_coef.begin();
		//auto ftr = b_coef.back();

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
    
    std::vector<double> matrix_pde_case1::backward_coefficient(const double& temp)
    {
		double dx_2 = pow(dx, 2.0);
		std::vector<double> c_coef(space_dim - 2);

		for (auto it = c_coef.begin(); it != c_coef.end(); ++it)
		{
			*it = temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		return c_coef;
    }
    
	void matrix_pde_case1::set_initial_conditions() //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change this
	{
		
		option_payoff.resize(space_dim - 1, 0.0);
		//filling  and initial condition in result vector
		
		option_payoff = pde->init_cond(S_values);
	}

	std::vector<double> matrix_pde_case1::boundary_increment(const double& t)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double> boundary_effect;
		boundary_effect.resize(space_dim -1, 0.0);
		//auto it1 = boundary_effect.begin(); //check this with prof
		//auto it2 = boundary_effect.back();
		double left_b = pde->boundary_left(t);
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
    
	dauphine::matrix matrix_pde_case1::transition_matrix(const double& temp)
	{
		std::vector<double> a = forward_coefficient(temp);
		std::vector<double> b = present_coefficient(temp);
		std::vector<double> c = backward_coefficient(temp);
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

	void matrix_pde_case1::Crout_Algo_Resolution()
	{
		double temp_lhs = theta * dt;
		double temp_rhs = - (1- theta) * dt;

		//setting initial conditions
		set_initial_conditions();


		//creation of the left and right transition matrices
		dauphine::matrix M_lhs = transition_matrix(temp_lhs);
		dauphine::matrix M_rhs = transition_matrix(temp_rhs);

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

	std::vector<double> matrix_pde_case1::LU_compute(dauphine::matrix& L, dauphine::matrix& U, const std::vector<double>& b)
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

	std::vector<double> matrix_pde_case1::get_option_payoff()
	{
		return option_payoff;
	}

	std::vector<double> matrix_pde_case1::get_price_curve()
	{
		//test if resolution of the PDE has been done
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		return new_result;
	}

	double matrix_pde_case1::get_price(const double& S)
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

	double matrix_pde_case1::compute_delta(const double& S)
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

	std::vector<double> matrix_pde_case1::compute_delta()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> delta(space_dim - 1);
		for (size_t i = 0; i < space_dim - 1;++i) 
		{
			delta.push_back(compute_delta(S_values[i]));
		}
		
		
		return delta;
	}

	double matrix_pde_case1::compute_gamma(const double& S)
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

	std::vector<double> matrix_pde_case1::compute_gamma()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> gamma(space_dim - 1);
		for (size_t i = 0; i < space_dim - 1; ++i)
		{
			gamma.push_back(compute_gamma(S_values[i]));
		}
		return gamma;
	}

	double matrix_pde_case1::compute_theta(const double& S)
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
		return (old_result[i] - new_result[i]) / dt * 365; //theta per 1 day change
	}

	std::vector<double> matrix_pde_case1::compute_theta()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		std::vector<double> theta(space_dim - 1);
		for (size_t i = 0; i < space_dim - 1; ++i)
		{
			theta.push_back(compute_theta(S_values[i]));
		}
		return theta;
	}
	
	std::vector<double> matrix_pde_case1::compute_vega()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}

		BS_PDE* vega_pde = pde->vega_pde();
		matrix_pde_case1* vega_solve = new Solve::matrix_pde_case1(vega_pde, theta, space_dim, time_dim, S0, maturity);
		vega_solve->Crout_Algo_Resolution();
		std::vector<double> vega_price = vega_solve->get_price_curve();
		vega.resize(space_dim - 1);
		for (size_t i = 0; i < space_dim - 1; ++i) {
			vega.push_back(vega_price[i] - new_result[i]);
		}
		return vega;
	}

	double matrix_pde_case1::compute_vega(const double& S)
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



