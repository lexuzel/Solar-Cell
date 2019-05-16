#pragma once
#include <vector>
#include "constants.h"

class Builder {
protected:
	std::vector<double> n_table;
	std::vector<double> alfa_table;

	double width, Ndop, D_diff, mu, tau, S;
	double step_x, step_l;

	/* ���������� ������� ����������� */
	/* ��������� ��� �������� ������� ������������ */
	virtual std::vector<double> construct_coef_matrix();

	/* ���������� ������ֲ� � ���ֲ. */
	/* ������� ����'���� calc_abs_coef(�); integrate_abs_coef(0, �) */
	virtual double calc_Generation(double x);

	/* ������������ ���Բֲ���� ���������� �� � */
	/* ������� ����'���� calc_abs_coef(x) */
	virtual double integrate_abs_coef(double start, double end);

	/* ���������� ���Բֲ���� ���������� � ���ֲ */
	/* ³��������� ������������ �� �������� ����� � ������ K_l
	   ������� ����'���� get_alfa(lambda) ��
	   get_solar_distribution(lambda) */
	virtual double calc_abs_coef(double x);

	virtual void fill_alfa_table();

	/* ���������� dndx � ���ֲ */
	virtual double calc_dndx(double x);

	virtual double get_alfa(double lambda, double x);

	virtual double get_Band_gap(double x);
	/* ���������� ����������Ҳ ���. ���� � ���ֲ */
	virtual double calc_Eq(double x);

	double get_solar_distribution(double lambda);
	bool double_equal(double a, double b);

public:
	/* в������ �����������Ҳ. */
	/* ����'���� ������� ��������� ������ � ������ K_x.
	   �� �����: n(x), dndx(x);
	   ���������� ��� ���������� ������� ����������� ��
	   ��������� construct_coef_matrix();
	   ������� ����'���� calculate_Generation(x)*/
	virtual void integrate_continuity_eq();

	/* �������� ��Ͳ � ���� */
	virtual void write_to_file(const char* filename);
};

class PZoneBuilder : public Builder {
protected:
	virtual std::vector<double> construct_coef_matrix() override final;
	virtual double calc_Generation(double x) override final;
	virtual double integrate_abs_coef(double start, double end) override final;
	virtual double calc_abs_coef(double x) override final;
	virtual void fill_alfa_table() override final;
	virtual double calc_dndx(double x) override final;
	virtual double get_alfa(double lambda, double x) override final;

	virtual double get_Band_gap(double x);
	virtual double calc_Eq(double x);
public:	
	PZoneBuilder();
	virtual void integrate_continuity_eq() override final;
	virtual void write_to_file(const char* filename) override final;
};

class NZoneBuilder : public Builder {
protected:
	const double delta_alfa;

	virtual std::vector<double> construct_coef_matrix() override final;
	virtual double calc_Generation(double x) override final;
	virtual double integrate_abs_coef(double start, double end) override final;
	virtual double calc_abs_coef(double x) override final;
	virtual void fill_alfa_table() override final;
	virtual double calc_dndx(double x) override final;
	virtual double get_alfa(double lambda, double x) override final;

	virtual double get_Band_gap(double x);
	virtual double calc_Eq(double x);
public:
	NZoneBuilder(double d_alfa);
	virtual void integrate_continuity_eq();
	virtual void write_to_file(const char* filename) override final;
};

class GradedPZoneBuilder : public PZoneBuilder {
protected:
	virtual double calc_Eq(double x) override final;
	virtual double get_Band_gap(double x) override final;
public:
	GradedPZoneBuilder();
};

class GradedNZoneBuilder : public NZoneBuilder {
protected:
	virtual double calc_Eq(double x) override final;
	virtual double get_Band_gap(double x) override final;
public:
	GradedNZoneBuilder(double d_alfa);
};

class SolarCell {
private:
	Builder builder;
public:
	void setBuilder(const Builder &b);
	void execute();
	void writeToFile(const char* address);
};