

#include "pde_solver.hpp"

namespace Solve
{
	Exo_Solver::Exo_Solver(Exo_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity)
		:BS_Solver(_pde, _theta, _space_dim, _time_dim, _S0, _maturity), pde(_pde), theta(_theta), space_dim(_space_dim),
		time_dim(_time_dim), S0(_S0), maturity(_maturity), dt(0.0), dx(0.0), x_min(0.0), x_max(0.0),resolved(false) {
		calculate_parameters();
	}
	
	void Exo_Solver::calculate_parameters() {
		double stdv = pde->standard_dev();

		x_max = log(S0) + 5 * stdv;
		x_min = log(S0) - 5 * stdv;
		dx = 10 * stdv / static_cast<double>(space_dim);
		dt = maturity / static_cast<double>(time_dim);
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

	void Exo_Solver::set_initial_conditions()
	{

		option_payoff.resize(space_dim - 1, 0.0);
		//filling  and initial condition in result vector

		option_payoff = pde->init_cond(S_values);
	}

	std::vector<double> Exo_Solver::forward_coefficient(const double& temp, const size_t& i)
	{
		Exo_Solver::calculate_parameters();
		double dx_2 = pow(dx, 2.0);
		std::vector<double>a_coef(space_dim - 2);
		std::vector<double>&& diff_coeff = pde->diff_coeff();
		std::vector<double>&& conv_coeff = pde->conv_coeff(i);

		for (size_t i = 0; i != space_dim - 2; ++i)
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
		Exo_Solver::calculate_parameters();

		double temp_lhs = theta * dt;
		double temp_rhs = -(1 - theta) * dt;

		//setting initial conditions
		Exo_Solver::set_initial_conditions();

		dauphine::matrix M_lhs(space_dim - 1, space_dim - 1);
		dauphine::matrix M_rhs(space_dim - 1, space_dim - 1);
		dauphine::matrix L(space_dim - 1, space_dim - 1);
		dauphine::matrix U(space_dim - 1, space_dim - 1);
		std::vector<double> v(space_dim - 1);
		std::vector<double> tmp(space_dim - 1);

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

			for (std::size_t i = 1; i != time_dim ; ++i)
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
			tmp = boundary_increment(time_dim - i);
			v = M_rhs.produit_mat_vect(old_result);
			std::transform(tmp.begin(), tmp.end(), v.begin(), v.begin(), std::plus<double>());

			new_result = LU_compute(L, U, v);
		}
		resolved = true;
		return;
	}

	std::vector<double> Exo_Solver::get_price_curve()
	{
		//test if resolution of the PDE has been done
		if (resolved == 0) {
			Crout_Algo_Resolution();
		}
		return new_result;
	}

	double Exo_Solver::get_price(const double& S)
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

	double Exo_Solver::compute_delta(const double& S)
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

	std::vector<double> Exo_Solver::compute_delta()
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

	double Exo_Solver::compute_gamma(const double& S)
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
			ds_2 = (S_values[space_dim - 2] - S_values[space_dim - 3]) * (S_values[space_dim - 3] - S_values[space_dim - 4]) * (S_values[space_dim - 2] - S_values[space_dim - 4]);
			return 2 / ds_2 * ((new_result[space_dim - 2] - new_result[space_dim - 3]) * (S_values[space_dim - 3] - S_values[space_dim - 4]) - (new_result[space_dim - 3] - new_result[space_dim - 4]) * (S_values[space_dim - 2] - S_values[space_dim - 3]));
		}
		else {
			ds_2 = (S_values[i + 1] - S_values[i]) * (S_values[i] - S_values[i - 1]) * (S_values[i + 1] - S_values[i - 1]);
			return 2 / ds_2 * ((new_result[i + 1] - new_result[i]) * (S_values[i] - S_values[i - 1]) - (new_result[i] - new_result[i - 1]) * (S_values[i + 1] - S_values[i]));
		}
	}

	std::vector<double> Exo_Solver::compute_gamma()
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

	double Exo_Solver::compute_theta(const double& S)
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

	std::vector<double> Exo_Solver::compute_theta()
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

	std::vector<double> Exo_Solver::get_option_payoff()
	{
		return option_payoff;
	}







}



