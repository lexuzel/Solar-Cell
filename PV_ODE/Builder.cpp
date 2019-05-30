#include "pch.h"

double Builder::calc_abs_coef(double x){
	double lambda = 5.5E-7;
	return get_alfa(lambda, x);;
}

double Builder::get_alfa(double lambda, double x){
	double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };

	if (d_energy < (get_Band_gap(x) / q_e)) return alfa0 * exp(d_energy - (get_Band_gap(x) - q_e));
	else return alfa0 * sqrt(d_energy);
}

double Builder::get_Band_gap(double x){
	return Eg_0 * q_e;
}

double Builder::calc_Eq(double x){
	return 0.0;
}

double Builder::get_solar_distribution(double lambda) {
	int base_index{ 0 };
	for (int i = 0; i < irra_table_size; i++) {
		if (X[i] > lambda * 1.0E6) {
			base_index = i - 1;
			break;
		}
		if (double_equal(X[i], lambda * 1.0E6)) {
			return Y[i];
		}
	}
	double k{ (Y[base_index + 1] - Y[base_index]) / (X[base_index + 1] - X[base_index]) };
	return Y[base_index] + k * (lambda * 1.0E6 - X[base_index]);
}

bool Builder::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}

double Builder::get_Boundary_x1(){
	return 2.0E21;
}

double Builder::get_Boundary_x2(){
	return 0.0;
}
