#pragma once

/* Файл, в якому містяться константи програми */

/* ФІЗИЧНІ КОНСТАНТИ */

const double pi = 3.14159265359;
const double plank_h = 6.626E-34;		// Стала Планка		
const long speed_of_light = 299792458;	// Швидкість світла
const double q_e = 1.6021766208E-19;	// Заряд електрона (модуль)
const double m0 = 9.1093826E-31;
const double k_b = 1.38064852E-23;		// Стала Больцмана
const double T = 300.0;					// Температура в Кельвінах

// Ширина бази - абсорбуючої області [m]
const double width_n = 0.78E-6;			// d = 0.78 mkm (Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>)
const double width_p = 0.78E-6;

// Концентрація домішок [1/m3]
const double NA = 2.0E18;				// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>
const double ND = 2.0E18;

// Константи дифузії електронів і дірок [m2/(V*s)]
const double D_diff_n = 1.01E-3;		// Dn = 10.1 cm2/s (Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>)
const double D_diff_p = 1.01E-3;

// Рухливості електронів і дірок [m/s * m/V]
const double mu_n = 0.04;				// mu(n) = 400 cm2/(V*s) (Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>)
const double mu_p = 0.004;

// Час життя електронів і дірок 
const double tau_n = 1.0E-9;			// L(diff)n = sqrt(Dn * tau) = 3.16E-5 cm (Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>)
const double tau_p = 1.0E-9;

// Швидкість поверхневої рекомбінації
const double S_n = 1.0E3;					// S = 1.0E5-E6 (Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>)
const double S_p = 1.0E3;

// Ширина забороненої зони (крайні точки)
const double Eg_0 = 1.5 + 1.26;
const double Eg_d = 1.5;
const double Eg_d2 = 1.5;
const double convex = 0.0;

// Гранична умова на концентрацію при x = width (зона p-n переходу)
const double m_e = 0.064 * Eg_d * m0;
const double m_h = 0.35 * m0;

// Константа у коефіцієнті поглинання alfa0
//const double alfa0 = 1.0E6;				// Theoretical analysis of solar cells based on graded band-gap structures <G.Sassi>

/* КОНСТАНТИ ЧИСЛОВОГО МЕТОДУ */

const long K_x = 1000;					// Сітка дискретизація для х (координата в базі)
const long K_l = 400;					// Сітка дискретизація для lambda (довжина хвилі світла, що поглинається)

/* ПОХІДНІ КОНСТАНТИ */

const double step_x_n = width_n / K_x;		// Крок сітки для х
const double step_x_p = width_p / K_x;
const double step_l = 0.98E-6 / K_l;		// Крок сітки для lambda

/* ТАБЛИЦЯ СПЕКТРАЛЬНОГО РОЗПОДІЛУ РАДІАЦІЇ СОНЦЯ (АМ 1) */

const double X[]{ 0.0, 0.11, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43,
				0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69,
				0.71, 0.73, 0.75, 0.77, 0.79, 0.83, 0.88, 0.93, 0.98 };

const double Y[]{ 0.0, 0.001, 0.035, 0.055, 0.155, 0.36, 0.515, 0.77, 0.835, 0.91, 0.87, 1.335,
				1.35, 1.55, 1.55, 1.45, 1.4, 1.35, 1.4, 1.35, 1.35, 1.25, 1.2, 1.2, 1.1, 1.05,
				1.0, 0.95, 0.9, 0.85, 0.85, 0.76, 0.68, 0.6, 0.56 };

const int irra_table_size = 35;