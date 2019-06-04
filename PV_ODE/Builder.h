#pragma once
#include <vector>

class Builder {
protected:
	std::vector<double> n_table;
	std::vector<double> alfa_table;
	std::vector<double> integral_alfa_table;
	std::vector<double> generation_table;

	double width, Ndop, D_diff, mu, tau, S;
	double step_x, x0;

	/* ���������� ������� ����������� */
	/* ��������� ��� �������� ������� ������������ */
	virtual std::vector<double> construct_coef_matrix();

	/* в������ �����������Ҳ. */
	/* ����'���� ������� ��������� ������ � ������ K_x.
	   �� �����: n(x), dndx(x);
	   ���������� ��� ���������� ������� ����������� ��
	   ��������� construct_coef_matrix();
	   ������� ����'���� calculate_Generation(x)*/
	virtual void integrate_continuity_eq(double b1, double b2);

	/* ���������� ������ֲ� � ���ֲ. */
	/* ������� ����'���� calc_abs_coef(�); integrate_abs_coef(0, �) */
	double calc_Generation(double x);
	virtual double get_alfa(double x, double lambda);
	virtual double get_Band_gap(double x);

	/* ���������� ����������Ҳ ���. ���� � ���ֲ */
	virtual double calc_Eq(double x);
	virtual double calc_dEdx(double x);

	/* ���������� dndx � ���ֲ */
	double calc_dndx(double x);
	bool double_equal(double a, double b);
public:
	double photo_current;

	virtual void write_to_file(const char* filename) = 0;
	virtual void calc_photo_carriers() = 0;
	std::vector<double> calc_delta();
	void calc_photoCurrent();

};