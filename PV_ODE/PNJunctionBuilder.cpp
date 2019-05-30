#include "pch.h"

std::vector<double> PNJunctionBuilder::construct_coef_matrix()
{
	return std::vector<double>();
}

double PNJunctionBuilder::calc_Generation(double x)
{
	return 0.0;
}

void PNJunctionBuilder::fill_alfa_table()
{
}

double PNJunctionBuilder::calc_dndx(double x)
{
	return 0.0;
}

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

double PNJunctionBuilder::get_Band_gap(double x){
	auto energy = get_Energy_level(x);
	return energy.first - energy.second;
}

double PNJunctionBuilder::calc_Eq(double x){
	if (double_equal(x, x_1)) 
		return (get_Band_gap(x + step_x) - Eg_x1) / (2 * step_x);
	if (double_equal(x, x_2)) 
		return (Eg_x2 - get_Band_gap(x - step_x)) / (2 * step_x);

	return (get_Band_gap(x + step_x) - get_Band_gap(x - step_x)) / (2 * step_x);
}

PNJunctionBuilder::PNJunctionBuilder(double d_alfa) : delta_alfa(d_alfa){
	width = x_2 - x_1;
	step_x = step_x_pn;

	n_table.resize(K_x + 1);
	alfa_table.resize(2 * K_x + 1);
}

void PNJunctionBuilder::integrate_continuity_eq()
{
}

void PNJunctionBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);
	fout << "  x  \t  Conduction level  \t  Valence level\n\n";
	for (double x = x_1; x <= x_2; x += step_x_pn) {
		fout.scientific;
		fout.precision(4);
		fout.width(10);

		auto energy = get_Energy_level(x);
		fout << x << "\t" << energy.first / q_e << "\t" << energy.second / q_e << "\n";
	}
	fout.close();
}

double PNJunctionBuilder::integrate_abs_coef(double start, double end)
{
	return 0.0;
}

double PNJunctionBuilder::calc_delta()
{
	return 0.0;
}

double PNJunctionBuilder::get_majority_carriers(double x)
{
	return 0.0;
}
