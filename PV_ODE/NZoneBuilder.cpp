#include "pch.h"

NZoneBuilder::NZoneBuilder(std::vector<double> &d_alfa) {
	width = x_d - x_2;
	Ndop = ND;
	D_diff = D_diff_p;
	mu = mu_p;
	tau = tau_p;
	S = S_p;
	step_x = step_x_p;

	n_table.resize(K_x + 1); 
	alfa_table.resize(2001);
	std::copy(d_alfa.begin(), d_alfa.end(), integral_alfa_table.begin());
}

void NZoneBuilder::integrate_continuity_eq(){

	std::vector<double> coef_matrix = construct_coef_matrix();

		double alfa3{ 0 }, beta{ 0 };
		alfa3 = 2 * D_diff - mu * step_x * calc_Eq(x_2 + width - step_x);
		beta = S / D_diff - mu / D_diff * calc_Eq(x_2 + width - step_x);
		coef_matrix[0] += alfa3;
		coef_matrix[1] -= 2 * step_x * alfa3 * beta;
		coef_matrix[2] = 0;

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
		n_table[0] = get_Boundary_x2();
		if (n_table[0] != 0.0) n_table[0] = -n_table[0];
		for (long i = K_x - 2; i >= 0; i--) {
			n_table[(K_x - 2) - i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[(K_x - 2) - i]) / coef_matrix[i * 4 + 1];
		}
		n_table[K_x] = n_table[K_x - 2] - n_table[K_x - 1] * 2 * step_x * beta;
}

void NZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);
	fout << "  x  \t  p(x)  \t  dndx(x)  \t  (maj)n(x)  \t Alfa  \t  Conduction level  \t Valence level  \t  EQ Field\n\n";
		for (double x = x_2; x < x_2 + width || double_equal(x, x_2 + width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round((x - x_2) / step_x));
			fout << x << "\t" << n_table[index_n] << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << get_majority_carriers(x) << "\t" << get_alfa(x, 5.5e-7) << "\t";
			fout << get_Band_gap(x) / q_e - bias << "\t" << -bias << "\t" << calc_Eq(x) << "\n";
		}
	fout.close();
}

std::vector<double> NZoneBuilder::construct_coef_matrix(){
	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq(x_2 + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x / tau;
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq(x_2 + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 3] = -calc_Generation(x_2 + width - (i + 1) * step_x) * 2 * step_x*step_x;
	}
	return coef_matrix;
}

double NZoneBuilder::get_majority_carriers(double x){
	long double ni = 2 * pow(2 * pi * sqrt(m_e * m_h) * k_b * T, 1.5) / pow(plank_h, 3) * exp(-get_Band_gap(x) / (2 * k_b * T));

	int index_n = static_cast<int>((x - x_2) / step_x);
	long double result = Ndop * Ndop*exp(2 * -bias * q_e / k_b / T) / (n_table[index_n] + ni);
	return static_cast<double>(result);
}

double NZoneBuilder::calc_dndx(double x){
	int index = static_cast<int>(round((x - x_2) / step_x));
	assert(index >= 1 && index <= K_x - 1 && "In calc_dndx");
	return (n_table[index + 1] - n_table[index - 1]) / (2 * step_x);
}

double NZoneBuilder::get_Band_gap(double x){
	return Eg_x2 * q_e;
}

