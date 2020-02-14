#ifndef pde_cpp
#define pde_cpp

#include "pde.hpp"
#include <math.h>

BS_PDE::BS_PDE(VanillaOption* _option, const std::string& _left_boundary_type, const std::string& _right_boundary_type)
	: option(_option), left_boundary_type(_left_boundary_type), right_boundary_type(_right_boundary_type) {}

std::string BS_PDE::get_right_boundary_type() const {
	return right_boundary_type;
}

std::string BS_PDE::get_left_boundary_type() const {
	return left_boundary_type;
}

double BS_PDE::diff_coeff() const {
	double vol = option->sigma;
	return 0.5 * vol * vol;
}

double BS_PDE::conv_coeff() const {
	return (option->r) - diff_coeff();
}

double BS_PDE::zero_coeff() const {
	return (option->r);
}

double BS_PDE::boundary_left() const {
	return 0.00001;
}

double BS_PDE::boundary_right(const double& t, const double& x) const {
	double res = 1.0;

	if (left_boundary_type.compare("D") == 0) {
		res = (x - (option->K) * exp(-(option->r) * ((option->T) - t))); //careful to use exp(x) when calling this function
	}
	return res;
}

double BS_PDE::init_cond(const double& x) {
	return option->pay_off->operator()(x); //careful to use exp(x) when calling this function
}

std::vector<double> BS_PDE::init_cond(const std::vector<double>& X)
{
	size_t l = X.size();
	std::vector<double> res(l);
	//std::transform(X.begin(), X.end(), res.begin(), [this](double x)->{ return init_cond(x);});
	std::transform(X.begin(), X.end(), res.begin(),
		[this](double arg) { return BS_PDE::init_cond(arg); });
	return res;
}

double BS_PDE::standard_dev() {
	double vol = option->sigma;
	double maturity = option->T;
	return vol * sqrt(maturity);
}

BS_PDE* BS_PDE::vega_pde(const double& d)
{
	VanillaOption* vega_op = option->Option_vega(d);
	BS_PDE* vega_pde = new BS_PDE(vega_op, left_boundary_type, right_boundary_type);
	return vega_pde;
}

Exo_PDE::Exo_PDE(ExoticOption* _option, const std::string& _left_boundary_type , const std::string& _right_boundary_type)
	: BS_PDE(_option,_left_boundary_type, _right_boundary_type), option(_option), left_boundary_type(_left_boundary_type), right_boundary_type(_right_boundary_type) {}

std::vector<double> Exo_PDE::diff_coeff() 
{
	std::vector<double>&& vol = option->get_vol_TS();
	std::transform(vol.begin(), vol.end(), vol.begin(), 
		[](double d) { return 0.5 * d * d; });
	return vol;
}

std::vector<double> Exo_PDE::conv_coeff(const size_t & t)
{
        // Implementation: why using rvalue reference here?
        // std::move would be better and would avoid potential
        // dangling reference if you decide to change the return type
        // of get_yield_curve
	std::vector<double>&& R = option->get_yield_curve();
	std::vector<double>&& diff = diff_coeff();
	double r = R[t];
	std::transform(diff.begin(), diff.end(), diff.begin(),
		[r](double d) { return r - d;});

	return diff;
}

std::vector<double> Exo_PDE::zero_coeff()
{
        // return option->get_yield_curve() would allow
        // return value optimization from the compiler
	std::vector<double>&& R = option->get_yield_curve();
	return R;
}

double Exo_PDE::standard_dev()
{
	std::vector<double>&& vol = option->get_vol_TS();
	double maturity = option->T;
	return vol[(vol.size() + 1 )/2] * sqrt(maturity); //Volatility for S0
}

double Exo_PDE::boundary_right(const size_t& i, const double& dt, const double& x) const
{
	double res = 1.0;
	std::vector<double>&& R = option->get_yield_curve();
	if (left_boundary_type.compare("D") == 0) {
		res = (x - (option->K) * exp(-(R[i]) * ((option->T) - i * dt))); //careful to use exp(x) when calling this function
	}
	return res;
}



#endif // !1
