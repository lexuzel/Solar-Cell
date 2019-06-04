#pragma once
#include <vector>

class Builder {
protected:
	std::vector<double> n_table, eq_table;
	std::vector<double> alfa_table, integral_alfa_table;
	std::vector<double> generation_table;

	double width, D_diff, mu, tau;
	double step_x, x0;

	/* Ïîáóäóâàòè ìàòğèöş êîåô³ö³ºíò³â */
	/* Íåîáõ³äíî äëÿ âèğ³øåííÿ ğ³âíÿííÿ íåïåğåğâíîñò³ */
	virtual std::vector<double> construct_coef_matrix();

	/* Ğ²ÂÍßÍÍß ÍÅÏÅĞÅĞÂÍÎÑÒ². */
	virtual void integrate_continuity_eq(double b1, double b2);
	virtual void calc_equilibrium(double b1, double b2);

	/* ÎÁ×ÈÑËÅÍÍß ÃÅÍÅĞÀÖ²¯ Â ÒÎ×Ö². */
	double calc_Generation(double x);
	double get_alfa(double x, double lambda);

	/* ÏÎÁÓÄÎÂÀ ÅÍÅĞÃÅÒÈ×ÍÎ¯ Ä²ÀÃĞÀÌÈ. */
	virtual std::pair<double, double> get_Energy_level(double x) = 0;
	double get_Band_gap(double x);
	virtual double calc_Eq(double x) = 0;
	double calc_dEdx(double x);

	double calc_dndx(double x);
	bool double_equal(double a, double b);
public:
	double photo_current;

	virtual void write_to_file(const char* filename) = 0;
	virtual void calc_photo_carriers() = 0;
	void calc_photoCurrent();

};