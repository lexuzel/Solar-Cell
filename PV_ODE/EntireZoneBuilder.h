#pragma once
#include "Builder.h"

class EntireZoneBuilder : public Builder {
	std::vector<double> electrons;
protected:
	virtual double calc_Eq(double x) override final;
	virtual std::pair<double, double> get_Energy_level(double x) override final;
	void set_ntype_constants();
	void set_ptype_constants();
public:
	EntireZoneBuilder();
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers() override final;
};