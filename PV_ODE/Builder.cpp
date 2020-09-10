#include "pch.h"

std::vector<double> Builder::integrate_continuity_eq(double xStart, double xEnd, double b1, double b2, bool equilibrium, int grid, double step, double t){

	if (grid == K_x)
		step = width / K_x;

	std::vector<double> coef_matrix = construct_coef_matrix(xStart, xEnd, equilibrium, grid, step, t);
	coef_matrix[3] -= b1 * coef_matrix[2];

	double alfa3{ 0 }, beta{ 0 };
	if (!equilibrium) {
		alfa3 = 2 * D_diff - mu * step * get_Eq(step);
		beta = S / D_diff - mu / D_diff * get_Eq(step);
		coef_matrix[0] += alfa3;
		coef_matrix[1] -= 2 * step * alfa3 * beta;
	}

	for (long i = 0; i < grid - 2; i++) {
		coef_matrix[(i + 1) * 4 + 1] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4] / coef_matrix[i * 4 + 1];
		coef_matrix[(i + 1) * 4 + 3] -= coef_matrix[(i + 1) * 4 + 2] * coef_matrix[i * 4 + 3] / coef_matrix[i * 4 + 1];
	}
	std::vector<double> table(grid + 1);
	table[grid] = b2;
	for (long i = grid - 2; i >= 0; i--) {
		table[i + 1] = (coef_matrix[i * 4 + 3] - coef_matrix[i * 4] * table[i + 2]) / coef_matrix[i * 4 + 1];
	}
	if (!equilibrium)
	{
		double r1 = table[1] * 2 * step * beta;
		table[0] = (table[2] - r1 > 0) ? table[2] - r1 : 0.0;
	}
	else table[0] = b1;
	return table;
}

std::vector<double> Builder::construct_coef_matrix(double xStart, double xEnd, bool equilibrium, int grid, double step, double t){

	std::vector<double> coef_matrix((grid - 1) * 4);
	if (equilibrium) 
	{
		for (long i = 0; i < grid - 1; i++) {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step * get_Eq((i + 1)  * step + xStart);
			coef_matrix[i * 4 + 1] = -4 * D_diff + 2 * step * step * mu * get_dEdx((i + 1) * step + xStart);
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step * get_Eq((i + 1) * step + xStart);
			coef_matrix[i * 4 + 3] = 0.0;
		}
	}
	else 
	{
		if (generation_table[0] == 0) {
			for (double x = xStart; x < xEnd || double_equal(x, xEnd); x += step) {
				int index = static_cast<int>(round((x - x0) / step));
				generation_table[index] = calc_Generation(x);
			}
		}

		static const int index_wp = static_cast<int>(round(w_p / step_x));
		for (long i = 0; i < grid - 1; i++) {
			coef_matrix[i * 4 + 0] = 2 * D_diff + mu * step * get_Eq((i + 1) * step + xStart);
			coef_matrix[i * 4 + 1] = -4 * D_diff - 2 * step * step * (t / tau - mu * get_dEdx((i + 1) * step + xStart));
			coef_matrix[i * 4 + 2] = 2 * D_diff - mu * step * (get_Eq((i + 1) * step + xStart));
			if (D_diff == D_diff_n && i < index_wp ||
				D_diff == D_diff_p && i > index_wp)
				coef_matrix[i * 4 + 3] = -generation_table[i + 1] * 2 * step * step;
			else
				coef_matrix[i * 4 + 3] = 0.0;

		}
	}
	return coef_matrix;
}

double Builder::calc_Generation(double x)
{
	std::string spectrumFile = "SolarSpectrum-AM15-NREL.txt";

	auto getSolarIrradiation = [](std::string file)
	{
		std::vector< std::pair<double, double> > distributionVector;
		std::ifstream fin(file);
		std::string lambda, irra;
		while (!fin.eof())
		{
			fin >> lambda >> irra;
			double l{ std::stod(lambda) };
			double ir{ std::stod(irra) };
			distributionVector.push_back({ l*1.0e-9, ir });
		}
		fin.close();
		return distributionVector;
	};

	static std::vector< std::pair<double, double> > distributionVector;
	if (distributionVector.empty())
	{
		distributionVector = getSolarIrradiation(spectrumFile);
		alfa_table.resize(distributionVector.size());
	}

	int index_x = static_cast<int>(round(x / step_x));
	if (integral_alfa_table[index_x].empty())
		integral_alfa_table[index_x].resize(distributionVector.size() - 1);

	alfa_table[0] = get_alfa(x, distributionVector[0].first);
	alfa_table[1] = get_alfa(x, distributionVector[1].first);
	if (index_x != 0) 
		integral_alfa_table[index_x].at(0) = integral_alfa_table[index_x - 1].at(0) + (alfa_table[1] + alfa_table[0]) / 2 * step_x;
	else
		integral_alfa_table[index_x].at(0) = (alfa_table[1] + alfa_table[0]) / 2 * step_x;

	double N0{ distributionVector[0].second * distributionVector[0].first / plank_h / speed_of_light };
	double f1{ N0 * alfa_table[1] * exp(-integral_alfa_table[index_x].at(0)) };
	double f2;

	double integral{ 0 };
	for (int i = 1; i < distributionVector.size() - 1; i++) {
		alfa_table[i+1] = get_alfa(x, distributionVector[i+1].first);
		if (index_x != 0) 
			integral_alfa_table[index_x].at(i) = integral_alfa_table[index_x - 1].at(i) + (alfa_table[i+1] + alfa_table[i]) / 2 * step_x;
		else
			integral_alfa_table[index_x].at(i) = (alfa_table[i + 1] + alfa_table[i]) / 2 * step_x;

		N0 = distributionVector[i].second * distributionVector[i].first / plank_h / speed_of_light;
		f2 = N0 * alfa_table[i+1] * exp(-integral_alfa_table[index_x].at(i));
		integral += (f1 + f2) / 2;
		f1 = f2;
	}

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

double Builder::get_derivative(std::vector<double> &table, int index, double step) 
{
	if (index == 0)
		return (table[index + 1] - table[index]) / (step);
	if (index == K_x)
		return (table[index] - table[index - 1]) / (step);

	return (table[index + 1] - table[index - 1]) / (2 * step);
}

double Builder::get_Band_gap(double x){
	auto energy = get_Energy_level(x);
	return energy.first - energy.second;
}

double Builder::get_dEdx(double x){
	return (get_Eq(x + step_x) - get_Eq(x - step_x)) / (2 * step_x);
}

bool Builder::double_equal(double a, double b) {
	double absEpsilon{ 1.0E-10 };
	double relEpsilon{ 1.0E-5 };
	double diff = fabs(a - b);
	if (diff < absEpsilon) return true;
	return diff < (fabs(a) > fabs(b) ? fabs(a) : fabs(b)) * relEpsilon;
}

