#include "pch.h"

double Builder::calc_Generation(double x){
	std::ifstream fin("SolarSpectrum-AM15-NREL.txt");
	auto get_solar_distribution = [&fin]() {
		if (fin.eof()) std::runtime_error("END OF FILE");
		std::string lambda, irra;
		fin >> lambda >> irra;
		double l{ std::stod(lambda) };
		double ir{ std::stod(irra) };
		return std::make_pair(l*1.0e-9, ir);
	};

	std::pair<double, double> solar;
	try {
		solar = get_solar_distribution();
	}
	catch (std::runtime_error &e) {
		std::cout << e.what() << "\n";
	}
	auto alfa{ get_alfa(x, solar.first) };
	if (alfa_table[0] != 0) integral_alfa_table[0] += (alfa + alfa_table[0]) / 2 * step_x;
	alfa_table[0] = alfa;
	double N0{ solar.second * solar.first / plank_h / speed_of_light };
	double f1{ N0 * alfa * exp(-integral_alfa_table[0]) };

	double integral{ 0 };
	for (int i = 1; i < 2001; i++) {
		try {
			solar = get_solar_distribution();
		}
		catch (std::runtime_error &e) {
			std::cout << e.what() << "\n";
		}
		alfa = get_alfa(x, solar.first);
		if (alfa_table[i] != 0) integral_alfa_table[i] += (alfa + alfa_table[i]) / 2 * step_x;
		alfa_table[i] = alfa;
		N0 = solar.second * solar.first / plank_h / speed_of_light;
		double f2{ N0 * alfa * exp(-integral_alfa_table[i]) };
		integral += (f1 + f2) / 2;
		f1 = f2;
	}
	fin.close();
	return integral;
}

double Builder::calc_dndx(double x){
	int index = static_cast<int>(round((x-x0) / step_x));
	assert(index >= 1 && index <= K_x - 1 && "In calc_dndx");
	return (n_table[index + 1] - n_table[index - 1]) / (2 * step_x);
}

std::vector<double> Builder::construct_coef_matrix(){

	if (generation_table[0] == 0) {
	for (double x = x0; x < x0 + width || double_equal(x, x0 + width); x += step_x) {
		int index = static_cast<int>(round((x - x0) / step_x));
		assert(index >= 0 && index <= K_x && "In construct_coef_matrix");
		generation_table[index] = calc_Generation(x);
	}
	}


	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
		coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq((i + 1) * step_x + x0);
		coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x * ( 1 / tau - mu * calc_dEdx((i + 1) * step_x + x0));
		coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq((i + 1) * step_x + x0);
		coef_matrix[i * 4 + 3] = -generation_table[i + 1] * 2 * step_x*step_x;
	}
	return coef_matrix;
}

double Builder::get_alfa(double x, double lambda){
	double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };

	if (d_energy <= 0.0) {
		double eps{ 0.01 };
		return alfa0 * exp(d_energy / eps);
	}
	else return alfa0 * (sqrt(d_energy) + 1);
}

double Builder::get_Band_gap(double x){
	return Eg_x1 * q_e;
}

double Builder::calc_Eq(double x){
	return 0.0;
}

double Builder::calc_dEdx(double x){
	return 0.0;
}

bool Builder::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}

void Builder::integrate_continuity_eq(double b1, double b2){

	std::vector<double> coef_matrix = construct_coef_matrix();
	coef_matrix[3] -= b1 * coef_matrix[2];

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
	n_table[K_x] = b2;
	for (long i = K_x - 2; i >= 0; i--) {
		n_table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[i + 2]) / coef_matrix[i * 4 + 1];
	}
	n_table[0] = b1;
}

std::vector<double> Builder::calc_delta(){
	return integral_alfa_table;
}

void Builder::calc_photoCurrent(){
	double integral{ 0 };
	for (int i = 0; i < generation_table.size() - 1; i++) {
		integral += (generation_table[i] + generation_table[i + 1]) / 2 * step_x;
	}
	photo_current = q_e * integral;
}

