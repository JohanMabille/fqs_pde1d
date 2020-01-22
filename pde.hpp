#pragma once

#ifndef pde_hpp
#define pde_hpp

#include "Option.hpp"
#include <string>
#include <vector>

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

class BS_PDE : public basicPDE {
	private:
		VanillaOption* option;
		std::string right_boundary_type;
		std::string left_boundary_type;

	public :
		BS_PDE(VanillaOption* _option, const std::string& _left_boundary_type = "D", const std::string& _right_boundary_type = "D");

		std::string get_right_boundary_type() const;
		std::string get_left_boundary_type() const;


		double diff_coeff(double t = 0.0, double x = 0.0, double v = 0.0) const;
		double conv_coeff(double t = 0.0, double x = 0.0, double v = 0.0) const;
		double zero_coeff(double t = 0.0, double x = 0.0, double v = 0.0) const;
		double source_coeff(double t = 0.0, double x = 0.0, double v = 0.0) const;

		double boundary_left(double t = 0.0, double x = 0.0,  double v = 0.0) const;
		double boundary_right(double t = 0.0, double x = 0.0, double v = 0.0) const;
		
		double init_cond(double x);
		std::vector<double> init_cond(std::vector<double> X);
		double standard_dev() ;

		BS_PDE* vega_pde();
		
		/*
		double diff_coeff2(double t, double x, double v) const = 0;
		double conv_coeff2(double t, double x, double v) const = 0;
		double boundary_left2(double t, double x, double v) const = 0;
		double boundary_right2(double t, double x, double v) const = 0;
		
		*/
		
};




#endif // !1

