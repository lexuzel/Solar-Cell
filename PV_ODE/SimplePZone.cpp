#include "pch.h"


SimplePZoneBuilder::SimplePZoneBuilder() : PZoneBuilder() {}

double SimplePZoneBuilder::calc_abs_coef(double x){
	double lambda = 5.5E-7;
	return get_alfa(lambda, x);
}

double SimplePZoneBuilder::get_alfa(double lambda, double x){
	double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };
	double alfa0{ 1.23E7 };

	if (d_energy < (get_Band_gap(x) / q_e)) return alfa0 * exp(d_energy - (get_Band_gap(x) - q_e));
	else return alfa0 * sqrt(d_energy);
}

/* -----------------------GradedSimplePZoneBuilder--------------------*/

GradedSimplePZoneBuilder::GradedSimplePZoneBuilder() : SimplePZoneBuilder() {}

double GradedSimplePZoneBuilder::calc_Eq(double x){
	double result = 2 * x * convex + (Eg_d - Eg_0) * q_e / width - convex * width;
	return result / q_e;
}

double GradedSimplePZoneBuilder::get_Band_gap(double x){
	double result = convex * x*x + ((Eg_d - Eg_0) * q_e / width - convex * width) * x + Eg_0 * q_e;
	return result;
}


