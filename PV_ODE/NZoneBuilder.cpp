#include "pch.h"

NZoneBuilder::NZoneBuilder(double d_alfa) : delta_alfa(d_alfa) {
	width = width_p;
	Ndop = ND;
	D_diff = D_diff_p;
	mu = mu_p;
	tau = tau_p;
	S = S_p;
	step_x = step_x_p;

	n_table.resize(K_x + 1); 
	alfa_table.resize(2 * K_x + 1);
}

void NZoneBuilder::integrate_continuity_eq(){
	fill_alfa_table();
	std::vector<double> coef_matrix = construct_coef_matrix();

		double alfa3{ 0 }, beta{ 0 };
		alfa3 = 2 * D_diff - mu * step_x * calc_Eq(width_n + width - step_x);
		beta = S / D_diff - mu / D_diff * calc_Eq(width_n + width - step_x);
		coef_matrix[0] += alfa3;
		coef_matrix[1] -= 2 * step_x * alfa3 * beta;
		coef_matrix[2] = 0;

	for (long i = 0; i < K_x - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
		n_table[0] = 0.0;
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
	fout << "  x  \t  n(x)  \t  dndx(x)  \t  Alfa  \t  Band gap  \t  EQ Field\n\n";
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
	fout.close();
}

std::vector<double> NZoneBuilder::construct_coef_matrix(){
	std::vector<double> coef_matrix((K_x - 1) * 4);
	for (long i = 0; i < K_x - 1; i++) {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step_x * calc_Eq(width_n + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step_x * step_x / tau;
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step_x * calc_Eq(width_n + width - (i + 1) * step_x);
			coef_matrix[i * 4 + 3] = -calc_Generation(width_n + width - (i + 1) * step_x) * 2 * step_x*step_x;
	}
	return coef_matrix;
}

double NZoneBuilder::calc_Generation(double x){
	int index = static_cast<int>(round(2 * (x - width_n) / step_x));
	assert(index >= 0 && index <= (2 * K_x) && "In calc_Generation");

	double alfa, integral_alfa;

	if (alfa_table[index] != 0.0) alfa = alfa_table[index];
	else alfa = calc_abs_coef(x);

	integral_alfa = integrate_abs_coef(width_n, x) + delta_alfa;

	double G0 = 1.0E15;
	return G0 * alfa * exp(-integral_alfa);			//!!!! ÒÐÅÁÀ ÂÑÒÀÂÈÒÈ G(0)
}

double NZoneBuilder::integrate_abs_coef(double start, double end){
	double x{ start };
	double f1{ calc_abs_coef(x) };
	double f2, f3;
	double integral{ 0 };
	while (x < end - step_x || double_equal(x, end - step_x)) {
		f2 = calc_abs_coef(x + step_x / 2);
		f3 = calc_abs_coef(x + step_x);
		integral += step_x * (f1 + 4 * f2 + f3) / 6;
		x += step_x;
		f1 = f3;
	}
	return integral;
}

double NZoneBuilder::calc_abs_coef(double x){
	int index = static_cast<int>(round(2 * (x - width_n) / step_x));
	assert(index >= 0 && index <= (2 * K_x) && "In calc_abs_coef");
	if (alfa_table[index] != 0.0) return alfa_table[index];

	double lambda{ step_l };
	double lambda_max{ plank_h * speed_of_light / get_Band_gap(x) };

	double norm = 0.0;
	for (double l = lambda; l <= lambda_max; l += step_l) {
		norm += get_solar_distribution(l);
	}

	double f1{ get_solar_distribution(lambda)*get_alfa(lambda, x) };
	double f2, f3;
	double integral{ 0.0 };
	while (lambda < lambda_max - step_l) {
		f2 = get_solar_distribution(lambda + step_l / 2)*get_alfa(lambda + step_l / 2, x);
		f3 = get_solar_distribution(lambda + step_l)*get_alfa(lambda + step_l, x);
		integral += (f1 + 4 * f2 + f3) / 6;
		lambda += step_l;
		f1 = f3;
	}
	return integral / norm;
}

void NZoneBuilder::fill_alfa_table(){
	int index{ 0 };
		for (double x = width_n; x < width_n + width || double_equal(x, width_n + width); x += step_x / 2) {
			alfa_table[index] = calc_abs_coef(x);
			index++;
		}
}

double NZoneBuilder::calc_dndx(double x){
	int index = static_cast<int>(round((x - width_n) / step_x));
	assert(index >= 1 && index <= K_x - 1 && "In calc_dndx");
	return (n_table[index + 1] - n_table[index - 1]) / (2 * step_x);
}

double NZoneBuilder::get_alfa(double lambda, double x){
	double d_energy{ (plank_h * speed_of_light / lambda - get_Band_gap(x)) / q_e };
	double alfa0{ 1.23E7 };
	return alfa0 * sqrt(d_energy);
}

double NZoneBuilder::get_Band_gap(double x){
	return Eg_d;
}

double NZoneBuilder::calc_Eq(double x){
	return 0.0;
}



/* -----------------------GradedNZoneBuilder--------------------*/
GradedNZoneBuilder::GradedNZoneBuilder(double d_alfa) : NZoneBuilder(d_alfa) {}

double GradedNZoneBuilder::calc_Eq(double x){
	double result;
		x -= width_n;
		result = 2 * x * convex + (Eg_d2 - Eg_d) * q_e / width - convex * width;
	return -result / q_e;
}

double GradedNZoneBuilder::get_Band_gap(double x){
	double result;
		x -= width_n;
		result = convex * x*x + ((Eg_d2 - Eg_d) * q_e / width - convex * width) * x + Eg_d * q_e;
	return result;
}