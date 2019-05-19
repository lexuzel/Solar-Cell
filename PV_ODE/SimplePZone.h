#pragma once
#include "PZoneBuilder.h"

class SimplePZoneBuilder : public PZoneBuilder {
protected:
	virtual double calc_abs_coef(double x) override final;
	virtual double get_alfa(double lambda, double x) override final;
public:
	SimplePZoneBuilder();
};

class GradedSimplePZoneBuilder : public SimplePZoneBuilder {
protected:
	virtual double calc_Eq(double x) override final;
	virtual double get_Band_gap(double x) override final;
public:
	GradedSimplePZoneBuilder();
};