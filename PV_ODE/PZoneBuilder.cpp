#include "pch.h"

PZoneBuilder::PZoneBuilder() {
	n_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);

	width = x_1;
	Ndop = NA;
	D_diff = D_diff_n;
	mu = mu_n;
	tau = tau_n;
	S = S_n;
	step_x = step_x_n;
}

std::vector<double> PZoneBuilder::construct_coef_matrix(){
	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
		coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq((i + 1) * step_x);
		coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x / tau;
		coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq((i + 1) * step_x);
		coef_matrix[i * 4 + 3] = -calc_Generation(step_x * (i + 1)) * 2 * step_x*step_x;
	}
	return coef_matrix;
}

//double PZoneBuilder::calc_Generation(double x) {
//	int index = static_cast<int>(round(2 * (x) / step_x));
//	assert(index >= 0 && index <= (2 * K_x) && "In calc_Generation");
//
//	double alfa, integral_alfa;
//
//	if (alfa_table[index] != 0.0) alfa = alfa_table[index];
//	else alfa = calc_abs_coef(x);
//
//	integral_alfa = integrate_abs_coef(0, x);
//
//	return G0 * alfa * exp(-integral_alfa);
//}

//double PZoneBuilder::integrate_abs_coef(double start, double end) {
//	double x{ start };
//	double f1{ calc_abs_coef(x) };
//	double f2, f3;
//	double integral{ 0 };
//	while (x < end - step_x || double_equal(x, end - step_x)) {
//		f2 = calc_abs_coef(x + step_x / 2);
//		f3 = calc_abs_coef(x + step_x);
//		integral += step_x * (f1 + 4 * f2 + f3) / 6;
//		x += step_x;
//		f1 = f3;
//	}
//	return integral;
//}

double PZoneBuilder::get_majority_carriers(double x){
	int index_n = static_cast<int>(round(x / step_x));
	long double result = Ndop + n_table[index_n];
	return static_cast<double>(result);
}

//void PZoneBuilder::fill_alfa_table(){
//	int index{ 0 };
//		for (double x = 0.0; x < width || double_equal(x, width); x += step_x / 2) {
//			alfa_table[index] = calc_abs_coef(x);
//			index++;
//		}
//}

double PZoneBuilder::calc_dndx(double x){
	int index = static_cast<int>(round(x / step_x));
	assert(index >= 1 && index <= K_x - 1 && "In calc_dndx");
	return (n_table[index + 1] - n_table[index - 1]) / (2 * step_x);
}

void PZoneBuilder::integrate_continuity_eq(){

	std::vector<double> coef_matrix = construct_coef_matrix();

	double alfa3{ 0 }, beta{ 0 };
	alfa3 = 2 * D_diff - mu * step_x * calc_Eq(step_x);
	beta = S / D_diff - mu / D_diff * calc_Eq(step_x);
	coef_matrix[0] += alfa3;
	coef_matrix[1] -= 2 * step_x * alfa3 * beta;
	coef_matrix[2] = 0;

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
	n_table[K_x] = get_Boundary_x1();
	for (long i = K_x - 2; i >= 0; i--) {
		n_table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[i + 2]) / coef_matrix[i * 4 + 1];
	}
	n_table[0] = n_table[2] - n_table[1] * 2 * step_x * beta;

	for (auto &a : alfa_table) a = 0;
	for (auto &a : integral_alfa_table) a = 0;
	double integral{ 0 };
	for (double x = 0.0; x < width - step_x || double_equal(x, width - step_x); x+=step_x) {
		integral += (calc_Generation(x) + calc_Generation(x + step_x)) / 2 * step_x;
	}
	std::cout << "Photo-Current = " << integral*q_e << "\n";
}

void PZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);
	for (auto &a : alfa_table) a = 0;
	for (auto &a : integral_alfa_table) a = 0;

	fout << "  x  \t  n(x)  \t  dndx(x)  \t (maj)p(x) \t  Alfa  \t  Conduction level  \t  Valence level  \t EQ Field   \t Generation\n\n";
		for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round(x / step_x));
			fout << x << "\t" << n_table[index_n] << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << get_majority_carriers(x) << "\t" << get_alfa(x, 5.5e-7) << "\t";
			fout << get_Band_gap(x) / q_e << "\t" << 0.0 << "\t" << calc_Eq(x) << "\t" <<calc_Generation(x) << "\n";
		}
	fout.close();
}

/* -----------------------GradedPZoneBuilder--------------------*/
GradedPZoneBuilder::GradedPZoneBuilder() : PZoneBuilder() {}

double GradedPZoneBuilder::calc_Eq(double x) {
	double result = (Eg_x1 - Eg_0) / width;
	return result;
}

double GradedPZoneBuilder::get_Band_gap(double x) {
	double result = (Eg_x1 - Eg_0) * q_e / width * x + Eg_0 * q_e;
	return result;
}

