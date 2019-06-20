#include "pch.h"

EntireZoneBuilder::EntireZoneBuilder(double width_p, double width_n) : 
					electrons(K_x + 1), 
					holes(K_x + 1),
					eq_electrons(K_x + 1),
					eq_holes(K_x + 1) {
	x0 = 0.0;
	w_p = width_p;
	width = width_p + width_n;
	step_x = width / K_x;

	generation_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);
}

std::pair<double, double> EntireZoneBuilder::get_Energy_level(double x){

	double energy_of_conduction_zone, energy_of_valence_zone;
	auto calc_zone_quadratic = [&x, this](const double &Ex1, const double &Ex2) {
		double a{ 0 }, b{ 0 }, c{ 0 };
		double dE_x1 = (Eg_x2 - Eg_0) / w_p;

		a = Ex2 - Ex1 + dE_x1 * ((w_p - x2) / 2 - w_p + x1);
		a /= (pow(w_p - x1, 2) - (w_p - x1) * (w_p - x2));
		if (x < w_p) {
			b = dE_x1 - 2 * a * x1;
			c = Ex1 - dE_x1 * x1 + a * x1*x1;
		}
		else {
			a = a * (w_p - x1) / (w_p - x2) + (dE_x1 / 2 / (w_p - x2));
			b = -2 * a * x2;
			c = Ex2 + a * x2*x2;
		}
		return (a * x*x + b * x + c) * q_e;
	};

	auto calc_zone = [&x, this](const double &Ex1, const double &Ex2) {
		return (Ex2 - Ex1) / (x2 - x1) * q_e * (x-x1) + Ex1 * q_e;
	};

	if (x < x1 || double_equal(x, x1)) {
		energy_of_conduction_zone = (Eg_x2 - Eg_0) * q_e / w_p * x + Eg_0 * q_e;
		energy_of_valence_zone = 0.0;
	}
	else if (x > x2) {
		energy_of_conduction_zone = (Eg_x2 - bias) * q_e;
		energy_of_valence_zone = -bias * q_e;
	}
	else {
		energy_of_conduction_zone = calc_zone_quadratic(get_Energy_level(x1).first / q_e, Eg_x2 - bias);
		energy_of_valence_zone = calc_zone_quadratic(0.0, -bias);
	}

	return std::make_pair(energy_of_conduction_zone, energy_of_valence_zone);
}

double EntireZoneBuilder::get_Eq(double x){
	if (D_diff == D_diff_n) {
		if (x < x1) return (Eg_x2 - Eg_0) / w_p;
		if (x > x2) return 0.0;
		return (get_Energy_level(x + step_x).first / q_e - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);
	}
	else {
		if (x < x1 || x > x2) return 0.0;
		return (get_Energy_level(x + step_x).second / q_e - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);
	}
}

void EntireZoneBuilder::calc_photo_carriers(double Eg0) {
	Eg_0 = Eg0;
	double width_OOZ = sqrt(2 * eps * eps0 / q_e * (bias) * (NA + ND) / NA / ND) / 10;
	double p_OOZ = (width_OOZ - eps * eps0 * (Eg_x2 - Eg_0) / (w_p *q_e * ND) / 10) / (NA/ND + 1);
	double n_OOZ = width_OOZ - p_OOZ;

	x1 = w_p - p_OOZ;
	x2 = w_p + n_OOZ;

	auto calc_current = [this] (std::vector<double> &table) {
		photo_current += get_photoCurrent();
		//recomb_current += get_recombCurrent(table);
		minor_current += get_minorCurrent(table);
	};
	auto set_ntype_constants = [this] {
		D_diff = D_diff_n;
		mu = mu_n;
		tau = tau_n;
	};
	auto set_ptype_constants = [this] {
		D_diff = D_diff_p;
		mu = -mu_p;
		tau = tau_p;
	};
	auto equilibrium = true;
	double boudary_electron = 4 * pow(2 * pi * sqrt(me_p * mh_p * m0 * m0) * k_b * T, 3) / pow(plank_h, 6) * exp(-Eg_x2 * q_e / (k_b * T));
	double boudary_hole = 4 * pow(2 * pi * sqrt(me_n * mh_n * m0 * m0) * k_b * T, 3) / pow(plank_h, 6) * exp(-Eg_x2 * q_e / (k_b * T));

	set_ntype_constants();
	electrons = integrate_continuity_eq(0.0, 0.0);
	calc_current(electrons);
	eq_electrons = integrate_continuity_eq(boudary_electron / NA, ND, equilibrium);
	equ_current += get_minorCurrent(eq_electrons);
	
	set_ptype_constants();
	holes = integrate_continuity_eq(0.0, 0.0);
//	calc_current(holes);
	minor_current += get_minorCurrent(holes);
	eq_holes = integrate_continuity_eq(NA, boudary_hole / ND, equilibrium);
	equ_current += get_minorCurrent(eq_holes);

	for (double &a : alfa_table) a = 0.0;
	for (double &a : integral_alfa_table) a = 0.0;
	for (double &a : generation_table) a = 0.0;
}

void EntireZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(15);

	auto get_field = [this](double x) {
		double electron_field, hole_field;
		if (x > 0.0 || x < width) {
			electron_field = (get_Energy_level(x + step_x).first / q_e - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);
			hole_field = (get_Energy_level(x + step_x).second / q_e - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);
		}
		else {
			electron_field = 0.0;
			hole_field = 0.0;
		}
		return std::make_pair(electron_field, hole_field);
	};

	fout << "x \t n(x) \t p(x) \t eq_n(x) \t eq_p(x) \t dndx \t dpdx \t Conduction level \t Valence level \t Electron-field \t Hole-field \t alfa \t Generation \n\n";
	for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
		fout.scientific;
		fout.precision(4);
		fout.width(15);
		int index = static_cast<int>(round((x - x0) / step_x));
		auto energy = get_Energy_level(x);
		fout << x << "\t" << electrons[index] << "\t" << holes[index] 
			      << "\t" << eq_electrons[index] << "\t" << eq_holes[index] << "\t";
		if (index == 0 || index == K_x) fout << 0.0 << "\t" << 0.0 << "\t";
		else fout << get_derivative(electrons, x) << "\t" << get_derivative(holes, x) << "\t";
		fout << energy.first / q_e << "\t" << energy.second / q_e << "\t" << get_field(x).first << "\t"
			 << get_field(x).second << "\t" << get_alfa(x, 750e-9) << "\t" << generation_table[index] << "\n";
	}
	fout.close();
}

double EntireZoneBuilder::get_photoCurrent(){
	double integral{ 0 };
	for (int i = 0; i < generation_table.size() - 1; i++) {
		integral += (generation_table[i] + generation_table[i + 1]) / 2 * step_x;
	}
	return q_e * integral;
}

double EntireZoneBuilder::get_recombCurrent(std::vector<double> &table){
	double integral{ 0 };
	for (int i = 0; i < table.size() - 1; i++) {
		integral += (table[i] + table[i + 1]) / 2 * step_x;
	}
	return q_e * integral / tau;
}

double EntireZoneBuilder::get_minorCurrent(std::vector<double>& table){
	double current = 1.0e10;
	double delta{ 2*200 * step_x };
	std::vector<double> dndx_table(2*200);
	int i = 0;
	if (D_diff == D_diff_n) {
		int index_x2 = static_cast<int>(round(x2 / step_x));
		for (double x = x2 - delta; x < x2; x += step_x) {
			dndx_table[i] = get_derivative(table, x);
			i++;
		}
		int index_dndx = std::max_element(dndx_table.begin(), dndx_table.end()) - dndx_table.begin();
		int index = index_dndx + index_x2 - 400;
		for (int i = -10; i < 11; i++) {
			double temp = abs(q_e * (mu * table[index + i] * get_Eq((index + i) * step_x) + D_diff * get_derivative(table, ((index + i) * step_x))));
			if (temp < current) current = temp;
		}
	}
	else {
		int index_x1 = static_cast<int>(round(x1 / step_x));
		for (double x = x1; x < x1 + delta; x += step_x) {
			dndx_table[i] = get_derivative(table, x);
			i++;
		}
		int index_dndx = std::min_element(dndx_table.begin(), dndx_table.end()) - dndx_table.begin();
		int index = index_dndx + index_x1;
		for (int i = -10; i < 11; i++) {
			double temp = abs(q_e * (mu * table[index + i] * get_Eq((index + i) * step_x) + D_diff * get_derivative(table, ((index + i) * step_x))));
			if (temp < current) current = temp;
		}
	}
	return current;
}

double EntireZoneBuilder::get_quantum_eff(std::vector<double>& table){
	std::vector<double> table_0tau = integrate_continuity_eq(0.0, 0.0, false, 0.0);
	double integral{ 0.0 }, integral_0tau{ 0.0 };
	for (int i = 0; i < table.size() - 1; i++) {
		integral += (table[i] + table[i + 1]) / 2 * step_x;
		integral_0tau += (table_0tau[i] + table_0tau[i + 1]) / 2 * step_x;
	}
	return integral / integral_0tau;
}





