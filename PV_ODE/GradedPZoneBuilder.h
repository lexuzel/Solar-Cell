#pragma once
#include "PZoneBuilder.h"

class GradedPZoneBuilder : public PZoneBuilder {
protected:
	virtual double calc_Eq(double x) override final;
	virtual double get_Band_gap(double x) override final;
public:
	GradedPZoneBuilder();
};