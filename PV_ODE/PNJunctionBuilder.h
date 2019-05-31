#pragma once
#include "Builder.h"

class PNJunctionBuilder : public Builder {
protected:
	virtual std::vector<double> construct_coef_matrix() override final;
	virtual double calc_dndx(double x) override final;
	virtual double get_Band_gap(double x) override final;
	virtual double calc_Eq(double x) override final;
	std::pair<double, double> get_Energy_level(double x);
public:
	PNJunctionBuilder(std::vector<double> d_alfa);
	virtual void integrate_continuity_eq() override final;
	virtual void write_to_file(const char* filename) override final;
	virtual double get_majority_carriers(double x) override final;
};