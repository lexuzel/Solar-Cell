#pragma once
#include <vector>

class Builder {
protected:
	std::vector<double> n_table;
	std::vector<double> alfa_table;

	double width, Ndop, D_diff, mu, tau, S;
	double step_x;

	/* ���������� ������� ����������� */
	/* ��������� ��� �������� ������� ������������ */
	virtual std::vector<double> construct_coef_matrix() = 0;

	/* ���������� ������ֲ� � ���ֲ. */
	/* ������� ����'���� calc_abs_coef(�); integrate_abs_coef(0, �) */
	virtual double calc_Generation(double x) = 0;

	/* ���������� ���Բֲ���� ���������� � ���ֲ */
	/* ³��������� ������������ �� �������� ����� � ������ K_l
	   ������� ����'���� get_alfa(lambda) ��
	   get_solar_distribution(lambda) */
	virtual double calc_abs_coef(double x);

	virtual void fill_alfa_table() = 0;

	/* ���������� dndx � ���ֲ */
	virtual double calc_dndx(double x) = 0;

	virtual double get_alfa(double lambda, double x);

	virtual double get_Band_gap(double x);
	/* ���������� ����������Ҳ ���. ���� � ���ֲ */
	virtual double calc_Eq(double x);

	double get_solar_distribution(double lambda);
	bool double_equal(double a, double b);

	double get_Boundary_x1();
	double get_Boundary_x2();

public:
	/* в������ �����������Ҳ. */
	/* ����'���� ������� ��������� ������ � ������ K_x.
	   �� �����: n(x), dndx(x);
	   ���������� ��� ���������� ������� ����������� ��
	   ��������� construct_coef_matrix();
	   ������� ����'���� calculate_Generation(x)*/
	virtual void integrate_continuity_eq() = 0;

	/* �������� ��Ͳ � ���� */
	virtual void write_to_file(const char* filename) = 0;

	/* ������������ ���Բֲ���� ���������� �� � */
/* ������� ����'���� calc_abs_coef(x) */
	virtual double integrate_abs_coef(double start, double end) = 0;

	virtual double calc_delta() = 0;
	virtual double get_majority_carriers(double x) = 0;
};