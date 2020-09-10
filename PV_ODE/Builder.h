#pragma once
#include <vector>
#include <mutex>

#include "constants.h"

class Builder {
protected:
	std::vector<double> alfa_table;
	std::vector< std::vector<double> > integral_alfa_table;
	std::vector<double> generation_table;

	double D_diff, mu, tau, S;
	double step_x, x0, x1, x2, width, w_p;
	double Eg_0;

	/* Ïîáóäóâàòè ìàòğèöş êîåô³ö³ºíò³â */
	/* Íåîáõ³äíî äëÿ âèğ³øåííÿ ğ³âíÿííÿ íåïåğåğâíîñò³ */
	virtual std::vector<double> construct_coef_matrix(double xStart, double xEnd, bool equilibrium, int grid, double step, double t);

	/* Ğ²ÂÍßÍÍß ÍÅÏÅĞÅĞÂÍÎÑÒ². */
	virtual std::vector<double> integrate_continuity_eq(double xStart, double xEnd, double b1, double b2, bool equilibrium = false, int grid = K_x, double step = 0.0, double t = 1.0);

	/* ÎÁ×ÈÑËÅÍÍß ÃÅÍÅĞÀÖ²¯ Â ÒÎ×Ö². */
	double calc_Generation(double x);
	double get_alfa(double x, double lambda);

	/* ÏÎÁÓÄÎÂÀ ÅÍÅĞÃÅÒÈ×ÍÎ¯ Ä²ÀÃĞÀÌÈ. */
	virtual std::pair<double, double> get_Energy_level(double x) = 0;
	double get_Band_gap(double x);
	virtual double get_Eq(double x) = 0;
	double get_dEdx(double x);

	double get_derivative(std::vector<double> &table, int index, double step);
	bool double_equal(double a, double b);

	virtual double get_photoCurrent() = 0;
	virtual double get_minorCurrent() = 0;
	virtual double get_recCurrent(std::vector<double>& table) = 0;
public:
	double photo_current = 0.0;
	double minor_current = 0.0;
	double rec_current = 0.0;
	double equ_current = 0.0;
	double maxEg0;

	virtual void write_to_file(const char* filename) = 0;
	virtual void calc_photo_carriers(double width_p, double width_n, double Eg0, double S) = 0;
};