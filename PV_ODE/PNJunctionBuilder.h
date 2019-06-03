#pragma once
#include "Builder.h"

class PNJunctionBuilder : public Builder {
	std::vector<double> electrons;
	std::vector<double> delta_alfa;
protected:
	virtual double get_Band_gap(double x) override final;
	virtual double calc_Eq(double x) override final;

	std::pair<double, double> get_Energy_level(double x);
	void set_ntype_constants();
	void set_ptype_constants();
public:
	PNJunctionBuilder(std::vector<double> d_alfa);
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers() override final;
};