#include "pch.h"

PZoneBuilder::PZoneBuilder() {
	n_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);
	generation_table.resize(K_x + 1);

	x0 = 0.0;
	width = x_1;
	Ndop = NA;
	D_diff = D_diff_n;
	mu = mu_n;
	tau = tau_n;
	step_x = step_x_n;
}

double PZoneBuilder::get_majority_carriers(double x){
	int index_n = static_cast<int>(round(x / step_x));
	long double result = Ndop + n_table[index_n];
	return static_cast<double>(result);
}

void PZoneBuilder::calc_photo_carriers(){
	integrate_continuity_eq(0.0, 0.0);
	calc_photoCurrent();
}

void PZoneBuilder::integrate_continuity_eq(double b1, double b2){

	std::vector<double> coef_matrix = construct_coef_matrix();

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

void PZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);

	fout << " x \t n(x) \t dndx(x) \t (maj)p(x) \t Alfa \t Conduction level \t Valence level \t Fermi level \t EQ Field \t Generation \n\n";
		for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round(x / step_x));
			fout << x << "\t" << n_table[index_n] << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << get_majority_carriers(x) << "\t" << get_alfa(x, 5.5e-7) << "\t";
			fout << get_Band_gap(x) / q_e << "\t" << 0.0 << "\t" << fermi_level << "\t" << calc_Eq(x) << "\t" << generation_table[index_n] << "\n";
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

