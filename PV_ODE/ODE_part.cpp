#include "pch.h"


Cell_Intrinsic::Cell_Intrinsic(double V, double d_alfa) :
	n_table(K_x + 1), alfa_table(2 * K_x + 1), voltage(V), delta_alfa(d_alfa) {

	if (d_alfa == 0.0) is_n = true;
	else is_n = false;
	if (is_n) {
		width = width_n;
		Ndop = NA;
		D_diff = D_diff_n;
		mu = mu_n;
		tau = tau_n;
		S = S_n;

		step_x = step_x_n;
	}
	else {
		width = width_p;
		Ndop = ND;
		D_diff = D_diff_p;
		mu = mu_p;
		tau = tau_p;
		S = S_p;

		step_x = step_x_p;
	}
}

double Cell_Intrinsic::calc_boundary() {
	long double ni = 2 * pow(2 * pi * sqrt(m_e * m_h) * k_b * T, 1.5) / pow(plank_h/(2*pi), 3) * exp(-Eg_d * q_e / (2 * k_b * T));
	long double r1 = ni * ni / Ndop * (exp(q_e * voltage / (k_b * T)) - 1);
	return static_cast<double>(r1);
}

void Cell_Intrinsic::integrate_continuity_eq() {
	fill_alfa_table();
	std::vector<double> coef_matrix = construct_coef_matrix();
	double beta = S / D_diff - mu / D_diff * calc_Eq(step_x);

	n_table[0] = 0.0;
	n_table[2] = coef_matrix[3] / (coef_matrix[0] + coef_matrix[1] / (2 * step_x * beta));
	n_table[1] = n_table[2] / (2 * step_x * beta);

	for (int i = 1; i < K_x - 1; i++) {
		n_table[i + 2] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4 + 1] * n_table[i + 1] - coef_matrix[i * 4 + 2] * n_table[i]) / coef_matrix[i * 4];
	}

//	double alfa3{ 0 }, beta{ 0 };
//	if (is_n) {
//		alfa3 = 2 * D_diff - mu * step_x * calc_Eq(step_x);
//		beta = S / D_diff - mu / D_diff * calc_Eq(step_x);
//		coef_matrix[0] += alfa3;
//		coef_matrix[1] -= 2 * step_x * alfa3 * beta;
//		coef_matrix[2] = 0;
//	}
//	else {
//		alfa3 = 2 * D_diff - mu * step_x * calc_Eq(width_n + width - step_x);
//		beta = S / D_diff - mu / D_diff * calc_Eq(width_n + width - step_x);
//		coef_matrix[0] += alfa3;
//		coef_matrix[1] -= 2 * step_x * alfa3 * beta;
//		coef_matrix[2] = 0;
//	}
//
//
//	for (long i = 0; i < K_x - 2; i++) {
//		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
//		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
//	}
//	if (is_n) {
//		n_table[K_x] = calc_boundary();
//		for (long i = K_x - 2; i >= 0; i--) {
//			n_table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[i + 2]) / coef_matrix[i * 4 + 1];
//		}
//		n_table[0] = n_table[2] - n_table[1] * 2 * step_x * beta;
//	}
//	else {
//		n_table[0] = calc_boundary();
////		if (n_table[0] != 0.0) n_table[0] = -n_table[0];
//		for (long i = K_x - 2; i >= 0; i--) {
//			n_table[(K_x - 2) - i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * n_table[(K_x - 2) - i]) / coef_matrix[i * 4 + 1];
//		}
//		n_table[K_x] = n_table[K_x - 2] - n_table[K_x - 1] * 2 * step_x * beta;
//	}



}

std::vector<double> Cell_Intrinsic::construct_coef_matrix() {
	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
		if (is_n) {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq((i + 1) * step_x);
			coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x / tau;
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq((i + 1) * step_x);
			coef_matrix[i * 4 + 3] = -calc_Generation(step_x * (i + 1)) * 2 * step_x*step_x;
		}
		else {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq(width_n + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x / tau;
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq(width_n + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 3] = -calc_Generation(width_n + width - (i + 1) * step_x) * 2 * step_x*step_x;
		}
	}

	return coef_matrix;
}
