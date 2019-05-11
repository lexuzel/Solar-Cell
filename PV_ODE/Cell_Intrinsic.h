#pragma once
#include <vector>

class Cell_Intrinsic {
private:
	double voltage;
	bool is_n;
	std::vector<double> n_table;
	std::vector<double> alfa_table;
	double delta_alfa;

	double width, Ndop, D_diff, mu, tau, S;
	double step_x;

	std::vector<double> construct_coef_matrix();
	double get_alfa(double lambda, double x);
	double get_solar_distribution(double lambda);
	double get_Band_gap(double x);
	void fill_alfa_table();

	bool double_equal(double a, double b);
public:
	Cell_Intrinsic(double V = 0.0, double d_alfa = 0.0);

	double calc_boundary();

	/* в������ �����������Ҳ. */
	/* ����'���� ������� ��������� ������ � ������ K_x.
	   �� �����: n(x), dndx(x);
	   ���������� ��� ���������� ������� ����������� ��
	   ��������� construct_coef_matrix();
	   ������� ����'���� calculate_Generation(x)*/
	void integrate_continuity_eq();

	/* ���������� ������ֲ� � ���ֲ. */
	/* ������� ����'���� calc_abs_coef(�); integrate_abs_coef(0, �) */
	double calc_Generation(double x);

	/* ������������ ���Բֲ���� ���������� �� � */
	/* ������� ����'���� calc_abs_coef(x) */
	double integrate_abs_coef(double start, double end);

	/* ���������� ���Բֲ���� ���������� � ���ֲ */
	/* ³��������� ������������ �� �������� ����� � ������ K_l
	   ������� ����'���� get_alfa(lambda) �� 
	   get_solar_distribution(lambda) */
	double calc_abs_coef(double x);
	
	/* ���������� ����������Ҳ ���. ���� � ���ֲ */
	double calc_Eq(double x);

	/* ���������� dndx � ���ֲ */
	double calc_dndx(double x);

	/* �������� ��Ͳ � ���� */
	void write_to_file(const char* filename);

};