#include "pch.h"

EntireZoneBuilder::EntireZoneBuilder(){
	x0 = 0.0;
	width = x_d;
	step_x = x_d / K_x;

	eq_table.resize(K_x + 1);
	n_table.resize(K_x + 1);
	electrons.resize(K_x + 1);
	generation_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);
}

std::pair<double, double> EntireZoneBuilder::get_Energy_level(double x){

	double energy_of_conduction_zone, energy_of_valence_zone;
	auto calc_zone = [&x](const double &Ex1, const double &Ex2) {
		double a{ 0 }, b{ 0 }, c{ 0 };
		double dE_x1 = (Eg_x1 - Eg_0) / x_1;

		a = Ex2 - Ex1 + dE_x1 * ((x_0 - x_2) / 2 - x_0 + x_1);
		a /= (pow(x_0 - x_1, 2) - (x_0 - x_1) * (x_0 - x_2));
		if (x < x_0) {
			b = dE_x1 - 2 * a * x_1;
			c = Ex1 - dE_x1 * x_1 + a * x_1*x_1;
		}
		else {
			a = a * (x_0 - x_1) / (x_0 - x_2) + (dE_x1 / 2 / (x_0 - x_2));
			b = -2 * a * x_2;
			c = Ex2 + a * x_2*x_2;
		}
		return (a * x*x + b * x + c) * q_e;
	};

	if (x < x_1) {
		energy_of_conduction_zone = (Eg_x1 - Eg_0) * q_e / x_1 * x + Eg_0 * q_e;
		energy_of_valence_zone = 0.0;
	}
	else if (x > x_2) {
		energy_of_conduction_zone = (Eg_x2 - bias) * q_e;
		energy_of_valence_zone = -bias * q_e;
	}
	else {
		energy_of_conduction_zone = calc_zone(Eg_x1, Eg_x2 - bias);
		energy_of_valence_zone = calc_zone(0.0, -bias);
	}

	return std::make_pair(energy_of_conduction_zone, energy_of_valence_zone);
}

double EntireZoneBuilder::calc_Eq(double x){
	if (D_diff == D_diff_n) {
		if (x < x_1) return (Eg_x1 - Eg_0) / x_1;
		if (x > x_2) return 0.0;

		if (double_equal(x, x_1))
			return (get_Energy_level(x + step_x).first / q_e - Eg_x1) / (2 * step_x);
		if (double_equal(x, x_2))
			return (-bias - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);

		return (get_Energy_level(x + step_x).first / q_e - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);
	}
	else {
		if (x < x_1 || x > x_2) return 0.0;
		if (double_equal(x, x_1))
			return (get_Energy_level(x + step_x).second / q_e - Eg_x1) / (2 * step_x);
		if (double_equal(x, x_2))
			return (-bias - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);

		return (get_Energy_level(x + step_x).second / q_e - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);
	}
}

void EntireZoneBuilder::set_ntype_constants(){
	D_diff = D_diff_n;
	mu = mu_n;
	tau = tau_n;
}

void EntireZoneBuilder::set_ptype_constants(){
	D_diff = D_diff_p;
	mu = -mu_p;
	tau = tau_p;
}

void EntireZoneBuilder::calc_photo_carriers() {
	set_ntype_constants();
	integrate_continuity_eq(0.0, 0.0);
	calc_photoCurrent();
	calc_equilibrium(0.0, ND);
	std::transform(n_table.begin(), n_table.end(), eq_table.begin(), electrons.begin(), [](double a, double b) {
		return a + b;
	});
	
	set_ptype_constants();
	integrate_continuity_eq(0.0, 0.0);
	calc_equilibrium(NA, 0.0);
	std::transform(n_table.begin(), n_table.end(), eq_table.begin(), n_table.begin(), [](double a, double b) {
		return a + b;
	});
}

void EntireZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);

	fout << "x \t n(x) \t p(x) \t Conduction level \t Valence level \t Band gap \t Generation \n\n";
	for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
		fout.scientific;
		fout.precision(4);
		fout.width(10);
		int index = static_cast<int>(round((x - x0) / step_x));
		auto energy = get_Energy_level(x);
		fout << x << "\t" << electrons[index] << "\t" << n_table[index] << "\t"
			<< energy.first / q_e << "\t" << energy.second / q_e << "\t" << get_Band_gap(x) /q_e << "\t" << generation_table[index] << "\n";
	}
	fout.close();
}