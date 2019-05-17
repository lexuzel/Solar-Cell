#include "pch.h"

double Builder::get_solar_distribution(double lambda) {
	int base_index{ 0 };
	for (int i = 0; i < irra_table_size; i++) {
		if (X[i] > lambda * 1.0E6) {
			base_index = i - 1;
			break;
		}
		if (double_equal(X[i], lambda * 1.0E6)) {
			return Y[i];
		}
	}
	double k{ (Y[base_index + 1] - Y[base_index]) / (X[base_index + 1] - X[base_index]) };
	return Y[base_index] + k * (lambda * 1.0E6 - X[base_index]);
}

bool Builder::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}