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

	/* РІВНЯННЯ НЕПЕРЕРВНОСТІ. */
	/* Розв'язок методом скінченних різниць з кроком K_x.
	   На виході: n(x), dndx(x);
	   Попередньо слід побудувати матрицю коефіцієнтів за
	   допомогою construct_coef_matrix();
	   Потребує розв'язку calculate_Generation(x)*/
	void integrate_continuity_eq();

	/* ОБЧИСЛЕННЯ ГЕНЕРАЦІЇ В ТОЧЦІ. */
	/* Потребує розв'язку calc_abs_coef(х); integrate_abs_coef(0, х) */
	double calc_Generation(double x);

	/* ІНТЕГРУВАННЯ КОЕФІЦІЄНТУ ПОГЛИНАННЯ ПО Х */
	/* Потребує розв'язку calc_abs_coef(x) */
	double integrate_abs_coef(double start, double end);

	/* ОБЧИСЛЕННЯ КОЕФІЦІЄНТУ ПОГЛИНАННЯ В ТОЧЦІ */
	/* Відбувається інтегрування по довжинам хвиль з кроком K_l
	   Потребує розв'язку get_alfa(lambda) та 
	   get_solar_distribution(lambda) */
	double calc_abs_coef(double x);
	
	/* ОБЧИСЛЕННЯ НАПРУЖЕНОСТІ ЕКВ. ПОЛЯ В ТОЧЦІ */
	double calc_Eq(double x);

	/* ОБЧИСЛЕННЯ dndx В ТОЧЦІ */
	double calc_dndx(double x);

	/* ЗАПИСАТИ ДАНІ В ФАЙЛ */
	void write_to_file(const char* filename);

};