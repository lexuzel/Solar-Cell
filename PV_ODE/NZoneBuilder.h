#pragma once
#include "Builder.h"

class NZoneBuilder : public Builder {
protected:
	const double delta_alfa;

	virtual std::vector<double> construct_coef_matrix() override final;
	virtual double calc_Generation(double x) override final;
	virtual double calc_abs_coef(double x) override;
	virtual void fill_alfa_table() override final;
	virtual double calc_dndx(double x) override final;
	virtual double get_alfa(double lambda, double x) override;

	virtual double get_Band_gap(double x) override;
	virtual double calc_Eq(double x) override;
public:
	NZoneBuilder(double d_alfa);
	virtual void integrate_continuity_eq() override;
	virtual void write_to_file(const char* filename) override final;
	virtual double integrate_abs_coef(double start, double end) override final;
};