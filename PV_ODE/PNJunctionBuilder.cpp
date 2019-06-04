#include "pch.h"


std::pair<double, double> PNJunctionBuilder::get_Energy_level(double x){

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

	double energy_of_conduction_zone = calc_zone(Eg_x1, Eg_x2 - bias);
	double energy_of_valence_zone = calc_zone(0.0, -bias);

	return std::make_pair(energy_of_conduction_zone, energy_of_valence_zone);
}

void PNJunctionBuilder::set_ntype_constants(){
	Ndop = NA;
	D_diff = D_diff_n;
	mu = mu_n;
	tau = tau_n;
}

void PNJunctionBuilder::set_ptype_constants(){
	Ndop = ND;
	D_diff = D_diff_p;
	mu = mu_p;
	tau = tau_p;
}

double PNJunctionBuilder::get_Band_gap(double x){
	auto energy = get_Energy_level(x);
	return energy.first - energy.second;
}

double PNJunctionBuilder::calc_Eq(double x){
	if (D_diff == D_diff_n) {
		if (double_equal(x, x_1))
			return (get_Energy_level(x + step_x).first / q_e - Eg_x1) / (2 * step_x);
		if (double_equal(x, x_2))
			return (-bias - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);

		return (get_Energy_level(x + step_x).first / q_e - get_Energy_level(x - step_x).first / q_e) / (2 * step_x);
	}
	else {
		if (double_equal(x, x_1))
			return -(get_Energy_level(x + step_x).second / q_e - Eg_x1) / (2 * step_x);
		if (double_equal(x, x_2))
			return -(-bias - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);

		return -(get_Energy_level(x + step_x).second / q_e - get_Energy_level(x - step_x).second / q_e) / (2 * step_x);
	}
}

double PNJunctionBuilder::calc_dEdx(double x){
	return (calc_Eq(x + step_x) - calc_Eq(x - step_x)) / (2 * step_x);
}

PNJunctionBuilder::PNJunctionBuilder(std::vector<double> d_alfa) : delta_alfa(d_alfa){
	x0 = x_1;
	width = x_2 - x_1;
	step_x = step_x_pn;

	n_table.resize(K_x + 1);
	electrons.resize(K_x + 1);
	generation_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);
	std::copy(d_alfa.begin(), d_alfa.end(), integral_alfa_table.begin());
}

void PNJunctionBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);

	fout << "x \t n(x) \t p(x) \t Conduction level \t Valence level \t Fermi level \t Generation \n\n";
	for (double x = x_1; x < x_2 || double_equal(x, x_2); x += step_x) {
		fout.scientific;
		fout.precision(4);
		fout.width(10);
		int index = static_cast<int>(round((x-x_1) / step_x));
		auto energy = get_Energy_level(x);
		fout << x << "\t" << electrons[index] << "\t" << n_table[index] << "\t" 
			<< energy.first / q_e << "\t" << energy.second / q_e << "\t" << fermi_level << "\t" << generation_table[index] << "\n";
	}
	fout.close();
}

void PNJunctionBuilder::calc_photo_carriers(){
	set_ntype_constants();
	integrate_continuity_eq(0.0, ND);
	std::copy(n_table.begin(), n_table.end(), electrons.begin());

	set_ptype_constants();
	integrate_continuity_eq(NA, 0.0);
	calc_photoCurrent();
}
