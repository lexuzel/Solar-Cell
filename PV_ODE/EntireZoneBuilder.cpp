#include "pch.h"
#include "EntireZoneBuilder.h"

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

void EntireZoneBuilder::integrate_continuity_eq(double b1, double b2){

	std::vector<double> coef_matrix = construct_coef_matrix();
	coef_matrix[3] -= b1 * coef_matrix[2];

	double alfa3{ 0 }, beta{ 0 };
	alfa3 = 2 * D_diff - mu * step_x * calc_Eq(step_x);
	beta = S / D_diff - mu / D_diff * calc_Eq(step_x);
	coef_matrix[0] += alfa3;
	coef_matrix[1] -= 2 * step_x * alfa3 * beta;

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
	n_table[K_x] = b2;
	for (long i = K_x - 2; i >= 0; i--) {
		n_table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[i + 2]) / coef_matrix[i * 4 + 1];
	}
	n_table[0] = n_table[2] - n_table[1] * 2 * step_x * beta;
}

double EntireZoneBuilder::get_Band_gap(double x){
	auto energy = get_Energy_level(x);
	return energy.first - energy.second;
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

double EntireZoneBuilder::calc_dEdx(double x){
	return (calc_Eq(x + step_x) - calc_Eq(x - step_x)) / (2 * step_x);
}

void EntireZoneBuilder::set_ntype_constants(){
	Ndop = NA;
	D_diff = D_diff_n;
	mu = mu_n;
	tau = tau_n;
}

void EntireZoneBuilder::set_ptype_constants(){
	Ndop = ND;
	D_diff = D_diff_p;
	mu = -mu_p;
	tau = tau_p;
}

void EntireZoneBuilder::calc_equilibrium(double b1, double b2){

	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
		coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq((i + 1) * step_x + x0);
		coef_matrix[i * 4 + 1] = -4 * D_diff + 2 * step_x * step_x * mu * calc_dEdx((i + 1) * step_x + x0);
		coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq((i + 1) * step_x + x0);
		coef_matrix[i * 4 + 3] = 0;
	}

	coef_matrix[3] -= b1 * coef_matrix[2];

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
	eq_table[K_x] = b2;
	for (long i = K_x - 2; i >= 0; i--) {
		eq_table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * eq_table[i + 2]) / coef_matrix[i * 4 + 1];
	}
	eq_table[0] = b1;
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