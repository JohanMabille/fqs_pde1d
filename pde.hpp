#ifndef pde_hpp
#define pde_hpp

#include "Option.hpp"
#include <string>
#include <vector>
#include "matrix.hpp"

//Two variable diffusion equation :

class basicPDE {
public:
	virtual double diff_coeff(double t, double x, double v) const = 0;
	virtual double conv_coeff(double t, double x, double v) const = 0;
	virtual double zero_coeff(double t, double x, double v) const = 0;
	virtual double source_coeff(double t, double x, double v) const = 0;
	/*
	virtual double diff_coeff2(double t, double x, double v) const = 0;
	virtual double conv_coeff2(double t, double x, double v) const = 0;
	virtual double boundary_left2(double t, double x, double v) const = 0;
	virtual double boundary_right2(double t, double x, double v) const = 0;
	*/
	virtual double boundary_left(double t, double x, double v) const = 0;
	virtual double boundary_right(double t, double x, double v) const = 0;
	//virtual double init_cond(double x) const = 0;
	//virtual double standard_dev() const = 0;

};

class BS_PDE {
private:
	VanillaOption* option;
	std::string right_boundary_type;
	std::string left_boundary_type;

public:
	BS_PDE();
	BS_PDE(VanillaOption* _option, const std::string& _left_boundary_type = "D", const std::string& _right_boundary_type = "D");

	std::string get_right_boundary_type() const;
	std::string get_left_boundary_type() const;


	virtual double diff_coeff() const;
	virtual double conv_coeff() const;
	virtual double zero_coeff() const;

	double boundary_left() const;
	virtual double boundary_right(const double& t = 0.0, const double& x = 0.0) const;

	double init_cond(const double& x);
	std::vector<double> init_cond(const std::vector<double>& X);
	double standard_dev();

	BS_PDE* vega_pde(const double& d);
};

class Exo_PDE : public BS_PDE {
private:
	ExoticOption* option;
	std::string right_boundary_type;
	std::string left_boundary_type;

public:
	Exo_PDE(ExoticOption* _option, const std::string& _left_boundary_type = "D", const std::string& _right_boundary_type = "D");

	std::vector<double> diff_coeff();
	std::vector<double> conv_coeff(const size_t& t);
	std::vector<double> zero_coeff();
	double standard_dev();

	//R vector will need to be inversed most likely !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double boundary_right(const size_t& i, const double& dt = 0.0, const double& x = 0.0) const;

	BS_PDE* vega_pde();


};


#endif // !1

