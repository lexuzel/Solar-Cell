#include "pch.h"


double Cell_Intrinsic::calc_Generation(double x) {
	int index;
	if (is_n) index = static_cast<int>(round(2 * x / step_x));
	else index = static_cast<int>(round(2 * (x - width_n) / step_x));
	assert(index >= 0 && index <= (2 * K_x) && "In calc_Generation");

	double alfa, integral_alfa;

	if (alfa_table[index] != 0.0) alfa = alfa_table[index];
	else alfa = calc_abs_coef(x);

	if (is_n) integral_alfa = integrate_abs_coef(0, x);
	else integral_alfa = integrate_abs_coef(width_n, x) + delta_alfa;

	double G0 = 1.0E15;

	return G0 * alfa * exp(-integral_alfa);			//!!!! ÒÐÅÁÀ ÂÑÒÀÂÈÒÈ G(0)
}

void Cell_Intrinsic::fill_alfa_table() {
	int index{ 0 };
	if (is_n) {
		for (double x = 0.0; x < width || double_equal(x, width); x += step_x / 2) {
			alfa_table[index] = calc_abs_coef(x);
			index++;
		}
	}
	else {
		for (double x = width_n; x < width_n + width || double_equal(x, width_n + width); x += step_x / 2) {
			alfa_table[index] = calc_abs_coef(x);
			index++;
		}
	}

}

double Cell_Intrinsic::integrate_abs_coef(double start, double end) {
	double x{ start };
	double f1{ calc_abs_coef(x) };
	double f2, f3;
	double integral{ 0 };
	while (x < end - step_x || double_equal(x, end - step_x)) {
		f2 = calc_abs_coef(x + step_x / 2);
		f3 = calc_abs_coef(x + step_x);
		integral += step_x * (f1 + 4 * f2 + f3) / 6;
		x += step_x;
		f1 = f3;
	}
	return integral;
}

double Cell_Intrinsic::calc_abs_coef(double x) {
	int index;
	if (is_n) index = static_cast<int>(round(2 * x / step_x));
	else index = static_cast<int>(round(2 * (x - width_n) / step_x));
	assert(index >= 0 && index <= (2 * K_x) && "In calc_abs_coef");
	if (alfa_table[index] != 0.0) return alfa_table[index];

	double lambda{ step_l };
	double lambda_max{ plank_h * speed_of_light / get_Band_gap(x) };

	double norm = 0.0;
	for (double l = lambda; l <= lambda_max; l += step_l) {
		norm += get_solar_distribution(l);
	}

	double f1{ get_solar_distribution(lambda)*get_alfa(lambda, x) };
	double f2, f3;
	double integral{ 0.0 };
	while (lambda < lambda_max - step_l) {
		f2 = get_solar_distribution(lambda + step_l / 2)*get_alfa(lambda + step_l / 2, x);
		f3 = get_solar_distribution(lambda + step_l)*get_alfa(lambda + step_l, x);
		integral += (f1 + 4 * f2 + f3) / 6;
		lambda += step_l;
		f1 = f3;
	}
	return integral / norm;
}

double Cell_Intrinsic::get_alfa(double lambda, double x)
{
	//double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };
	//double alfa0, eps;
	//if (d_energy < 0.004) {
	//	alfa0 = 2.63E6;
	//	eps = 9.86E-3;
	//}
	//else if (d_energy < 0.1) {
	//	alfa0 = 3.741E6;
	//	eps = 7.52E-2;
	//}
	//else {
	//	alfa0 = 1.23E7;
	//	return alfa0 * sqrt(d_energy);
	//}
	//return alfa0 * exp(d_energy / eps);

	//double alfa0 = 1.0E7;
	//double ph_energy = plank_h * speed_of_light / lambda;
	//if (ph_energy < (Eg_d * q_e)) return 0.0;
	//double Eg_eff = 3 * (Eg_0 - Eg_d) / 2 + Eg_d;
	//if (ph_energy < (Eg_eff * q_e)) return 2 * alfa0 / 3 * (ph_energy - Eg_d) / (Eg_0 - Eg_d);
	//else return alfa0;

	double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };
	double alfa0{ 1.23E7 };
	return alfa0 * sqrt(d_energy);
}
