#pragma once
#include <vector>

class Builder {
protected:
	std::vector<double> n_table;
	std::vector<double> alfa_table;

	double width, Ndop, D_diff, mu, tau, S;
	double step_x;

	/* Побудувати матрицю коефіцієнтів */
	/* Необхідно для вирішення рівняння неперервності */
	virtual std::vector<double> construct_coef_matrix() = 0;

	/* ОБЧИСЛЕННЯ ГЕНЕРАЦІЇ В ТОЧЦІ. */
	/* Потребує розв'язку calc_abs_coef(х); integrate_abs_coef(0, х) */
	virtual double calc_Generation(double x) = 0;

	/* ОБЧИСЛЕННЯ КОЕФІЦІЄНТУ ПОГЛИНАННЯ В ТОЧЦІ */
	/* Відбувається інтегрування по довжинам хвиль з кроком K_l
	   Потребує розв'язку get_alfa(lambda) та
	   get_solar_distribution(lambda) */
	virtual double calc_abs_coef(double x);

	virtual void fill_alfa_table() = 0;

	/* ОБЧИСЛЕННЯ dndx В ТОЧЦІ */
	virtual double calc_dndx(double x) = 0;

	virtual double get_alfa(double lambda, double x);

	virtual double get_Band_gap(double x);
	/* ОБЧИСЛЕННЯ НАПРУЖЕНОСТІ ЕКВ. ПОЛЯ В ТОЧЦІ */
	virtual double calc_Eq(double x);

	double get_solar_distribution(double lambda);
	bool double_equal(double a, double b);

	double get_Boundary_x1();
	double get_Boundary_x2();

public:
	/* РІВНЯННЯ НЕПЕРЕРВНОСТІ. */
	/* Розв'язок методом скінченних різниць з кроком K_x.
	   На виході: n(x), dndx(x);
	   Попередньо слід побудувати матрицю коефіцієнтів за
	   допомогою construct_coef_matrix();
	   Потребує розв'язку calculate_Generation(x)*/
	virtual void integrate_continuity_eq() = 0;

	/* ЗАПИСАТИ ДАНІ В ФАЙЛ */
	virtual void write_to_file(const char* filename) = 0;

	/* ІНТЕГРУВАННЯ КОЕФІЦІЄНТУ ПОГЛИНАННЯ ПО Х */
/* Потребує розв'язку calc_abs_coef(x) */
	virtual double integrate_abs_coef(double start, double end) = 0;

	virtual double calc_delta() = 0;
	virtual double get_majority_carriers(double x) = 0;
};