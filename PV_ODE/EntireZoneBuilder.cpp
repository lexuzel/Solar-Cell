#include "pch.h"

EntireZoneBuilder::EntireZoneBuilder()
{
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

	double biasV = k_b * T / q_e * log(NA*ND / intrCarrierConc(PZONE));

	if (x < x1 || double_equal(x, x1)) {
		energy_of_conduction_zone = (Eg_x2 - Eg_0) * q_e / w_p * x + Eg_0 * q_e;
		energy_of_valence_zone = 0.0;
	}
	else if (x > x2) {
		energy_of_conduction_zone = (Eg_x2 - biasV) * q_e;
		energy_of_valence_zone = -biasV * q_e;
	}
	else {
		energy_of_conduction_zone = calc_zone_quadratic(get_Energy_level(x1).first / q_e, Eg_x2 - biasV);
		energy_of_valence_zone = calc_zone_quadratic(0.0, -biasV);
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

void EntireZoneBuilder::getDepletionRegionCoord()
{
	double biasV = k_b * T / q_e * log(NA*ND / intrCarrierConc(PZONE));
	double width_OOZ = sqrt(2 * eps * eps0 / q_e * biasV * ( 1/NA + 1/ND ));

	maxEg0 = Eg_x2 + (NA / ND) / (eps * eps0 / (w_p *q_e * ND)) * width_OOZ;
	//Eg_0 = 0.995 * maxEg0;

	double p_OOZ = (width_OOZ - eps * eps0 * (Eg_x2 - Eg_0) / (w_p *q_e * ND)) / (NA / ND + 1);
	p_OOZ = p_OOZ > width_OOZ ? width_OOZ : p_OOZ;
	double n_OOZ = width_OOZ - p_OOZ;

	x1 = (w_p - p_OOZ) > 0.0 ? (w_p - p_OOZ) : 0.0;
	x2 = (w_p + n_OOZ) < width ? (w_p + n_OOZ) : width;
}

void EntireZoneBuilder::calc_photo_carriers(double width_p, double width_n, double Eg0, double S) {

	this->S = S;
	x0 = 0.0;
	w_p = width_p;
	width = width_p + width_n;
	step_x = width / K_x;

	generation_table.resize(K_x + 1);
	integral_alfa_table.resize(K_x + 1);

	Eg_0 = Eg0;
	getDepletionRegionCoord();

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
	double boudary_electron = intrCarrierConc(PZONE);
	double boudary_hole = intrCarrierConc(NZONE);
	eq_electrons.resize(K_x + 1);
	eq_holes.resize(K_x + 1);

	set_ntype_constants();
	electrons = integrate_continuity_eq(0.0, width, 0.0, 0.0);
	rec_current = x2 != w_p ? get_recCurrent(electrons) : 0.0;
	eq_electrons = integrate_continuity_eq(0.0, 0.0, boudary_electron / NA, ND, true);
	
	set_ptype_constants();
	holes = integrate_continuity_eq(0.0, width, 0.0, 0.0);
	rec_current += x2 != w_p ? get_recCurrent(holes) : 0.0;
	eq_holes = integrate_continuity_eq(0.0, 0.0, NA, boudary_hole / ND, true);

	//photo_current = get_photoCurrent();
	minor_current = x2 != w_p ? get_minorCurrent() : 0.0;
	equ_current = get_satCurrent();

	generation_table.clear();
	integral_alfa_table.clear();
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

	int index = 0;
	for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
		fout.scientific;
		fout.precision(4);
		fout.width(15);
		auto energy = get_Energy_level(x);
		fout << x << "\t" << electrons[index] << "\t" << holes[index] << "\t"
			 << eq_electrons[index] << "\t" << eq_holes[index] << "\t";
		
		if (index == 0 || index == K_x) fout << 0.0 << "\t" << 0.0 << "\t";
		else fout << get_derivative(electrons, index, step_x) << "\t" << get_derivative(holes, index, step_x) << "\t";
		fout << energy.first / q_e << "\t" << energy.second / q_e << "\t" << get_field(x).first << "\t"
			 << get_field(x).second << "\t" << get_alfa(x, 750e-9) << "\t" << generation_table[index] << "\n";

		index++;
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

double EntireZoneBuilder::getAverageField(double p1, double p2)
{
	int index_x1 = static_cast<int>(round(x1 / step_x));
	int index_x2 = static_cast<int>(round(x2 / step_x));

	static std::vector<double> fieldSumElectron(index_x2 + 1);
	static std::vector<double> fieldSumHole(K_x - index_x1 + 1);

	if (fieldSumElectron[index_x2] == 0)
	{
		double sum = 0.0;
		int index = index_x2;
		for (double x = x2; x > 0.0 || double_equal(x, 0.0); x -= step_x)
		{
			sum += -get_Eq(x);
			fieldSumElectron[index] = sum;
			index--;
		}
	}

	if (fieldSumHole[0] == 0)
	{
		double sum = 0.0;
		int index = 0;
		for (double x = x1; x < width || double_equal(x, width); x += step_x)
		{
			sum += -get_Eq(x);
			fieldSumHole[index] = sum;
			index++;
		}
	}

	int index_p1 = static_cast<int>(round(p1 / step_x));
	int index_p2 = static_cast<int>(round(p2 / step_x));
	double electronResult = fieldSumElectron[index_p1] / (index_x2 - index_p1);
	double holeResult = fieldSumHole[index_p2 - index_x1] / (index_p2 - index_x1);
	return fieldSumHole[0] == 0 ? electronResult : holeResult;
}

double EntireZoneBuilder::get_minorCurrent()
{
	int index_x1 = static_cast<int>(round(x1 / step_x));
	int index_x2 = static_cast<int>(round(x2 / step_x));

	double er1 = q_e * (electrons[index_x2]) * mu_n * -get_Eq(x2);
	std::vector<double> table;
	table.push_back(electrons[index_x2 - 1]);
	table.push_back(electrons[index_x2]);
	table.push_back(electrons[index_x2 + 1]);
	double er2 = q_e * D_diff_n * get_derivative(table, 1, step_x);
	double electronCurrent = er1 - er2;

	double hr1 = q_e * (holes[index_x1]) * mu_p * -get_Eq(x1);
	table.clear();
	if (index_x1 != 0)
		table.push_back(holes[index_x1 - 1]);
	table.push_back(holes[index_x1]);
	table.push_back(holes[index_x1 + 1]);
	double hr2 = q_e * D_diff_p * get_derivative(table, table.size() - 2, step_x);
	double holeCurrent = hr1 - hr2;

	return electronCurrent + holeCurrent;
}

double EntireZoneBuilder::intrCarrierConc(Zone zone)
{
	double carrierConc;
	switch (zone)
	{
	case PZONE:
		carrierConc = 4 * pow(2 * pi * sqrt(me_p * mh_p * m0 * m0) * k_b * T, 3) / pow(plank_h, 6) * exp(-Eg_x1 * q_e / (k_b * T));
		break;
	case NZONE:
		carrierConc = 4 * pow(2 * pi * sqrt(me_n * mh_n * m0 * m0) * k_b * T, 3) / pow(plank_h, 6) * exp(-Eg_x2 * q_e / (k_b * T));
	default:
		;
	}
	return carrierConc;
}

double EntireZoneBuilder::get_satCurrent()
{
	// Reference:
	// https://www.youtube.com/watch?v=kVgFZsnpoQI

	double np0 = intrCarrierConc(PZONE) / NA;
	double nn0 = intrCarrierConc(NZONE) / ND;

	double Jnsat = q_e * D_diff_n * np0 / sqrt(D_diff_n * tau_n);
	double Jpsat = q_e * D_diff_p * nn0 / sqrt(D_diff_p * tau_p);
	return Jnsat + Jpsat;
}

double EntireZoneBuilder::get_recCurrent(std::vector<double>& table){

	static int index_wp = static_cast<int>(round(w_p / step_x));

	int startIdx = (D_diff == D_diff_n) ? index_wp : 0;
	int endIdx = (D_diff == D_diff_n) ? K_x : index_wp;
	double integral{ 0.0 };
	for (int i = 0; i < K_x - 1; i++) 
	{
		integral += (table[i] + table[i + 1]) / 2 * step_x;
	}
	return integral * q_e / tau;
}





