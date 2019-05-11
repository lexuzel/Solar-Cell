#include "pch.h"



double Cell_Intrinsic::get_solar_distribution(double lambda) {
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

double Cell_Intrinsic::get_Band_gap(double x) {
	double result;
	if (is_n) result = convex * x*x + ((Eg_d - Eg_0) * q_e / width - convex * width) * x + Eg_0 * q_e;
	else {
		x -= width_n;
		result = convex * x*x + ((Eg_d2 - Eg_d) * q_e / width - convex * width) * x + Eg_d * q_e;
	}
	return result;
}


bool Cell_Intrinsic::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}

double Cell_Intrinsic::calc_Eq(double x) {
	double result;
	if (is_n) result = 2 * x * convex + (Eg_d - Eg_0) * q_e / width - convex * width;
	else {
		x -= width_n;
		result = 2 * x * convex + (Eg_d2 - Eg_d) * q_e / width - convex * width;
	}
	return -result / q_e;
}

double Cell_Intrinsic::calc_dndx(double x) {
	int index;
	if (is_n) index = static_cast<int>(round(x / step_x));
	else index = static_cast<int>(round((x - width_n) / step_x));
	assert(index >= 1 && index <= K_x - 1 && "In calc_dndx");
	return (n_table[index + 1] - n_table[index - 1]) / (2 * step_x);
}

void Cell_Intrinsic::write_to_file(const char* filename) {
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);
	fout << "  x  \t  n(x)  \t  dndx(x)  \t  Alfa  \t  Band gap  \t  EQ Field\n\n";
	if (is_n) {
		for (double x = 0.0; x < width || double_equal(x, width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round(x / step_x));
			fout << x << "\t" << n_table[index_n] << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << calc_abs_coef(x) << "\t";
			fout << get_Band_gap(x) / q_e << "\t" << calc_Eq(x) << "\n";
		}
	}
	else {
		for (double x = width_n; x < width_n + width || double_equal(x, width_n + width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round((x - width_n) / step_x));
			fout << x << "\t" << n_table[index_n] << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << calc_abs_coef(x) << "\t";
			fout << get_Band_gap(x) / q_e << "\t" << calc_Eq(x) << "\n";
		}
	}
	fout.close();
}
