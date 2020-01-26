

#include "pde_solver.hpp"

namespace Solve
{
	Exo_Solver::Exo_Solver(Exo_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity)
		:pde(_pde), theta(_theta), space_dim(_space_dim), time_dim(_time_dim), S0(_S0), maturity(_maturity) {
		calculate_parameters();
	}

	std::vector<double> Exo_Solver::forward_coefficient(const double& temp, const size_t& i)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double>a_coef(space_dim - 2);
		std::vector<double>&& diff_coeff = pde->diff_coeff();
		std::vector<double>&& conv_coeff = pde->conv_coeff(i);

		for (size_t i = 0; it != space_dim - 2; ++i)
		{
			a_coef[i] = temp * (conv_coeff[i] / (2 * dx) - diff_coeff[i] / dx_2);
		}
		return a_coef;
	}

	std::vector<double> Exo_Solver::present_coefficient(const double& temp, const size_t& i)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double>b_coef(space_dim - 1);
		std::vector<double>&& diff_coeff = pde->diff_coeff();
		std::vector<double>&& conv_coeff = pde->conv_coeff(i);
		double zero_coeff = pde->zero_coeff()[i];
		if (l.compare("N") == 0)
		{
			b_coef[0] = 1 + temp * (zero_coeff + 2 * diff_coeff[0] / dx_2) + temp * (conv_coeff[0] / (2 * dx) - diff_coeff[0] / dx_2);
		}
		else
		{
			b_coef[0] = 1 + temp * (zero_coeff + 2 * diff_coeff[0] / dx_2);
		}

		for (std::size_t i = 1; i < space_dim - 2; ++i)
		{
			b_coef[i] = 1 + temp * (zero_coeff + 2 * diff_coeff[i] / dx_2);
		}

		if (r.compare("N") == 0)
		{
			//present coeff + forward coeff because of neumann condition
			b_coef[space_dim - 2] = 1 + temp * (zero_coeff + 2 * diff_coeff[space_dim - 2] / dx_2) + temp * (conv_coeff[space_dim - 2] / (2 * dx) - diff_coeff[space_dim - 2] / dx_2);
		}
		else
		{
			b_coef[space_dim - 2] = 1 + temp * (zero_coeff + 2 * diff_coeff[space_dim - 2] / dx_2);
		}

		return b_coef;
	}

	std::vector<double> Exo_Solver::backward_coefficient(const double& temp, const size_t& i)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double> c_coef(space_dim - 2);
		std::vector<double>&& diff_coeff = pde->diff_coeff();
		std::vector<double>&& conv_coeff = pde->conv_coeff(i);

		for (size_t i = 0; i != space_dim - 2; ++i)
		{
			c_coef[i] = temp * (conv_coeff[i + 1] / (2 * dx) - diff_coeff[i + 1] / dx_2);
		}
		return c_coef;
	}

	std::vector<double> Exo_Solver::boundary_increment(const size_t& i)
	{
		double dx_2 = pow(dx, 2.0);
		std::vector<double> boundary_effect;
		boundary_effect.resize(space_dim - 1, 0.0);
		std::vector<double>&& diff_coeff = pde->diff_coeff();
		std::vector<double>&& conv_coeff = pde->conv_coeff(i);
		double left_b = pde->boundary_left();
		double right_b = pde->boundary_right(i, dt, exp(x_max));

		if (l.compare("D") == 0)
		{
			boundary_effect[0] = -left_b * dt * (conv_coeff[0] / (2 * dx) - diff_coeff[0] / dx_2); //Boudary times the backward coefficient.
		}
		else
		{
			boundary_effect[0] = left_b * dx * dt * (conv_coeff[0] / (2 * dx) - diff_coeff[0] / dx_2); //Boudary times the backward coefficient.
		}

		if (r.compare("D") == 0)
		{
			boundary_effect[space_dim - 2] = -right_b * dt * (-conv_coeff[space_dim - 2] / (2 * dx) - diff_coeff[space_dim - 2] / dx_2); //Boudary times the forward coefficient.
		}
		else
		{
			boundary_effect[space_dim - 2] = right_b * dx * dt * (-conv_coeff[space_dim - 2] / (2 * dx) - diff_coeff[space_dim - 2] / dx_2); //Boudary times the forward coefficient.
		}
		return boundary_effect;
	}

	void Exo_Solver::Crout_Algo_Resolution()
	{
		double temp_lhs = theta * dt;
		double temp_rhs = -(1 - theta) * dt;

		//setting initial conditions
		set_initial_conditions();

		dauphine::matrix&& M_lhs(space_dim - 1, space_dim - 1);
		dauphine::matrix&& M_rhs(space_dim - 1, space_dim - 1);
		dauphine::matrix L(space_dim - 1, space_dim - 1);
		dauphine::matrix U(space_dim - 1, space_dim - 1);
		std::vector<double> v(space_dim - 1);
		std::vector<double> tmp(space_dim - 1);
		double t;

		new_result = option_payoff;

		for (std::size_t i = 1; i != time_dim; ++i)
		{
			//creation of the left and right transition matrices
			M_lhs = transition_matrix(temp_lhs, i - 1);
			M_rhs = transition_matrix(temp_rhs, i - 1);

			// L - U decomposition of the left matrix for the Crout Algorithm 
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


			old_result = new_result;
			t = maturity - i * dt; // think about && rvalue lvalue;
			tmp = boundary_increment(t);
			v = M_rhs.produit_mat_vect(old_result);
			std::transform(tmp.begin(), tmp.end(), v.begin(), v.begin(), std::plus<double>());

			new_result = LU_compute(L, U, v);
		}
		resolved = true;
	}

	std::vector<double> BS_Solver::compute_vega()
	{
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		vega.resize(space_dim - 1);
		BS_PDE* vega_pde1 = pde->vega_pde(0.01);
		BS_PDE* vega_pde2 = pde->vega_pde(-0.01);
		BS_Solver* vega_solve1 = new Solve::BS_Solver(vega_pde1, theta, space_dim, time_dim, S0, maturity);
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
		if (vega.size() == 0) {
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



