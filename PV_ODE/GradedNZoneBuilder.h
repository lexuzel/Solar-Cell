#pragma once
#include "NZoneBuilder.h"

class GradedNZoneBuilder : public NZoneBuilder {
protected:
	virtual double calc_Eq(double x) override final;
	virtual double get_Band_gap(double x) override final;
public:
	GradedNZoneBuilder(double d_alfa);
};