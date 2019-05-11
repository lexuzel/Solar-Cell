
#include <iostream>
#include <vector>
#include <fstream>

const int grid = 1000;
const double width_1 = 0.78E-6;
const double width_2 = 0.5E-6;			//???
const double width_3 = 1.0E-5;

double correction = -width_1 - width_2 / 2;

double step_1 = width_1 / grid;
double step_2 = width_2 / grid;
double step_3 = width_3 / grid;

const double R = 0.0;				//???
const double tau = 1.0E-9;					
const double F0 = 1.0E21;			//???
const double gamma = 1.0E7;			//??? 1.0E7
const double diff = 1.01E-3;			// 1.0E-4
const double mu = 4.0E-3;
const double E1 = 1.615E6;

const double q = 1.6021766208E-19;
const double k = 1.38064852E-23;
const double T = 300.0;

long double theta = gamma * sqrt(diff*tau);
long double calc_pph = (1 - R)*(gamma * tau / (theta + 1))*F0;

long double Lp = sqrt(diff * tau);
long double d = mu * E1 * tau;
long double l = sqrt(d*d + 4 * diff*tau);
long double Z = exp(l * (width_1 + width_2 / 2) / (diff * tau));

long double L1 = (l + d) / 2;
long double L2 = (l - d) / 2;
long double eps2 = q * width_2 * E1 / (2 * k * T);
long double a = l / (Lp * (exp(l * width_1 / (diff * tau)) - 1));

long double A_r1 = (1 + a + Lp / L2 - q * E1 * Lp / (2 * k * T)) * sqrt(eps2) / (tan(sqrt(eps2)) * width_2);
long double A_r2 = 1 / L2 + a / Lp - q * E1 / (2 * k*T) - eps2 * Lp / (width_2 * width_2);
long double A = A_r1 + A_r2;

long double A1 = sqrt(eps2) * exp(-eps2 / 2) / (width_2 * sin(sqrt(eps2)));
long double K1 = A1 / A;
double K1_x1 = static_cast<double>(K1);

long double A2 = q*E1 / (2*k*T) - sqrt(eps2) / (width_2 * tan(sqrt(eps2))) + 1/L2 + a/Lp;
long double K2 = A2 / A;
double K2_x2 = static_cast<double>(K2);

std::vector<double> p_first(grid);
std::vector<double> p_second(grid);
std::vector<double> p_third(grid);

double psi(double x) {
	long double psi = Z * exp(x / L2) - exp(-x / L1);
	return static_cast<double>(psi);
}

double calc_K2(double x) {
	long double K2_r1 = (K2_x2 - 1 / (1 - theta)) * exp((width_1 + width_2 + correction - x) / Lp);
	long double K2_r2 = exp(gamma * (width_1 + width_2 + correction - x)) / (1 - theta);
	return static_cast<double>(K2_r1 + K2_r2);
}

double alfa(double x) {
	return sqrt(eps2) * (x / width_2);
}

double fi(double x) {
	return eps2 * pow((width_1 + width_2 + correction - x) / width_2, 2);
}

double fi_1(double x) {
	long double fi_r1 = sin(alfa(width_1 + width_2 + correction) - alfa(x)) / sin(sqrt(eps2));
	long double fi_r2 = exp((fi(width_1 + correction) - fi(x)) / 2);
	long double fi_1 = fi_r1 * fi_r2;
	return static_cast<double>(fi_1);
}

double fi_2(double x) {
	long double fi_r1 = -sin(alfa(width_1 + correction) - alfa(x)) / sin(sqrt(eps2));
	long double fi_r2 = exp((fi(width_1 + width_2 + correction) - fi(x)) / 2);
	long double fi_2 = fi_r1 * fi_r2;
	return static_cast<double>(fi_2);
}

void calc_zone1() {
	for (int i = 0; i < grid; i++) {
		p_first[i] = psi(i * step_1 + correction) / psi(width_1 + correction) * K1_x1 * calc_pph;
	}
}

void calc_zone2() {
	for (int i = 0; i < grid; i++) {
		p_second[i] = (K1_x1 * fi_1(i * step_2 + width_1 + correction) + K2_x2 * fi_2(i * step_2 + width_1 + correction)) * calc_pph;
	}
}

void calc_zone3() {
	for (int i = 0; i < grid; i++) {
		p_third[i] = calc_K2(i * step_3 + width_1 + width_2 + correction) * calc_pph;
	}
}

int main()
{
	calc_zone1();
	calc_zone2();
	calc_zone3();

	std::ofstream fout;
	fout.open("test.txt");
	fout << "x\t" << "p(x)\t" << "dp/dt\n";
	for (int i = 0; i < grid; i++) {
		fout << i * step_1 + correction << "\t" << p_first[i] << "\t";
		if (i == 0) fout << 0.0 << "\n";
		else if (i == grid - 1) fout << static_cast<double>((p_second[0] - p_first[grid - 2]) / (2 * step_1)) << "\n";
		else {
			double dndx = (p_first[i + 1] - p_first[i - 1]) / (2 * step_1);
			fout << dndx << "\n";
		}
	}
	fout << "\n";
	for (int i = 0; i < grid; i++) {
		fout << (i * step_2 + width_1) + correction << "\t" << p_second[i] << "\t";
		if (i == 0) fout << static_cast<double>((p_second[1] - p_first[grid - 1]) / (2 * step_2)) << "\n";
		else if (i == grid - 1) fout << static_cast<double>((p_third[0] - p_second[grid - 2]) / (2 * step_2)) << "\n";
		else {
			double dndx = (p_second[i + 1] - p_second[i - 1]) / (2 * step_2);
			fout << dndx << "\n";
		}
	}
	fout << "\n";
	for (int i = 0; i < grid; i++) {
		fout << (i * step_3 + width_1 + width_2) + correction << "\t" << p_third[i] << "\t";
		if (i == 0) fout << static_cast<double>((p_third[1] - p_second[grid - 1]) / (2 * step_3)) << "\n";
		else if (i == grid - 1) fout << 0.0 << "\n";
		else {
			double dndx = (p_third[i + 1] - p_third[i - 1]) / (2 * step_3);
			fout << dndx << "\n";
		}
	}
	fout.close();
	return 0;
}

