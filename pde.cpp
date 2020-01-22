#ifndef pde_cpp
#define pde_cpp

#include "pde.hpp"
#include <math.h>

BS_PDE::BS_PDE(VanillaOption* _option, const std::string& _left_boundary_type , const std::string& _right_boundary_type)
	: option(_option), left_boundary_type(_left_boundary_type), right_boundary_type(_right_boundary_type) {}

std::string BS_PDE::get_right_boundary_type() const {
	return right_boundary_type ;
}

std::string BS_PDE::get_left_boundary_type() const {
	return left_boundary_type ;
}

double BS_PDE::diff_coeff(double t, double x, double v) const {
	double vol = option->sigma;
	return 0.5 * vol * vol ;  
}

double BS_PDE::conv_coeff(double t, double x, double v) const {
	return (option->r) - diff_coeff(t, x, v) ;
}

double BS_PDE::zero_coeff(double t, double x, double v) const {
	return (option->r) ;
}

//maybe delete this
double BS_PDE::source_coeff(double t, double x, double v) const {
	return 0.0;
}

double BS_PDE::boundary_left(double t, double x, double v) const {
	return 0.0;
}

double BS_PDE::boundary_right(double t, double x, double v) const {
	double res = 1.0;

	if (left_boundary_type.compare("D") == 0) {
		res = (x - (option->K) * exp(-(option->r) * ((option->T) - t))); //careful to use exp(x) when calling this function
	}
	return res;
}

double BS_PDE::init_cond(double x)  {
	return option->pay_off->operator()(x); //careful to use exp(x) when calling this function
}

std::vector<double> BS_PDE::init_cond(std::vector<double> X)
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

BS_PDE* BS_PDE::vega_pde()
{
	VanillaOption* vega_op = option->Option_vega();
	BS_PDE* vega_pde = new BS_PDE(vega_op, left_boundary_type, right_boundary_type);
	return vega_pde;
}

#endif // !1
