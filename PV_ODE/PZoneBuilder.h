#pragma once
#include "Builder.h"

class PZoneBuilder : public Builder {
protected:
	virtual void integrate_continuity_eq(double b1, double b2) override final;
	double get_majority_carriers(double x);
public:
	PZoneBuilder();
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers() override final;
};