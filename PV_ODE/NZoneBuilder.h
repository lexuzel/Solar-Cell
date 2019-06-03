#pragma once
#include "Builder.h"

class NZoneBuilder : public Builder {
protected:
	std::vector<double> delta_alfa;

	virtual double get_Band_gap(double x) override final;
	double get_majority_carriers(double x);
public:
	NZoneBuilder(std::vector<double> &d_alfa);
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers() override final;
};