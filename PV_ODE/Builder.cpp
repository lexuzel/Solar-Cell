#include "pch.h"

//double Builder::calc_abs_coef(double x){
//	double lambda = 5.5E-7;
//	return get_alfa(x, lambda);
//}

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

//double Builder::get_solar_distribution(double lambda) {
//	int base_index{ 0 };
//	for (int i = 0; i < irra_table_size; i++) {
//		if (X[i] > lambda * 1.0E6) {
//			base_index = i - 1;
//			break;
//		}
//		if (double_equal(X[i], lambda * 1.0E6)) {
//			return Y[i];
//		}
//	}
//	double k{ (Y[base_index + 1] - Y[base_index]) / (X[base_index + 1] - X[base_index]) };
//	return Y[base_index] + k * (lambda * 1.0E6 - X[base_index]);
//}

bool Builder::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}

double Builder::get_Boundary_x1(){
	double correction = 0.5;
	long double ni = 2 * pow(2 * pi * sqrt(m_e * m_h) * k_b * T, 1.5) / pow(plank_h, 3) * exp(-get_Band_gap(x_1) * correction / (2 * k_b * T));
//	return static_cast<double>(ni*ni / NA);
	return 0.0;
}

double Builder::get_Boundary_x2(){
	return 0.0;
}

std::vector<double> Builder::calc_delta(){
	return integral_alfa_table;
}

