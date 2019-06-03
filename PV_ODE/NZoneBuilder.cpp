#include "pch.h"

NZoneBuilder::NZoneBuilder(std::vector<double> &d_alfa) : delta_alfa(d_alfa){
	x0 = x_2;
	width = x_d - x_2;
	Ndop = ND;
	D_diff = D_diff_p;
	mu = mu_p;
	tau = tau_p;
	step_x = step_x_p;

	n_table.resize(K_x + 1); 
	generation_table.resize(K_x + 1);
	alfa_table.resize(2001);
	integral_alfa_table.resize(2001);
	std::copy(d_alfa.begin(), d_alfa.end(), integral_alfa_table.begin());
}

void NZoneBuilder::write_to_file(const char * filename){
	std::ofstream fout;
	fout.open(filename);
	fout.width(10);

	fout << " x \t (maj)n(x) \t dndx(x) \t p(x) \t Alfa \t Conduction level \t Valence level \t Fermi level \t EQ Field \t Generation \n\n";
		for (double x = x_2; x < x_2 + width || double_equal(x, x_2 + width); x += step_x) {
			fout.scientific;
			fout.precision(4);
			fout.width(10);
			int index_n = static_cast<int>(round((x - x_2) / step_x));
			fout << x << "\t" << get_majority_carriers(x) << "\t";

			if (index_n > 0 && index_n < K_x) fout << calc_dndx(x) << "\t";
			else fout << 0.0 << "\t";

			fout << n_table[index_n] << "\t" << get_alfa(x, 5.5e-7) << "\t";
			fout << get_Band_gap(x) / q_e - bias << "\t" << -bias << "\t" << fermi_level << "\t" << calc_Eq(x) << "\t" << generation_table[index_n] << "\n";
		}
	fout.close();
}

double NZoneBuilder::get_majority_carriers(double x){
	int index_n = static_cast<int>(round((x - x_2) / step_x));
	long double result = Ndop + n_table[index_n];
	return static_cast<double>(result);
}

void NZoneBuilder::calc_photo_carriers(){
	integrate_continuity_eq(0.0, 0.0);
	calc_photoCurrent();
}

double NZoneBuilder::get_Band_gap(double x){
	return Eg_x2 * q_e;
}

