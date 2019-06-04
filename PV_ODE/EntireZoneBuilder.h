#pragma once
#include "Builder.h"

class EntireZoneBuilder : public Builder {
	std::vector<double> electrons;
	std::vector<double> eq_table;
protected:
	virtual void integrate_continuity_eq(double b1, double b2) override final;
	virtual double get_Band_gap(double x) override final;
	virtual double calc_Eq(double x) override final;
	virtual double calc_dEdx(double x) override final;
	std::pair<double, double> get_Energy_level(double x);
	void set_ntype_constants();
	void set_ptype_constants();
	void calc_equilibrium(double b1, double b2);
public:
	EntireZoneBuilder();
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers() override final;
};