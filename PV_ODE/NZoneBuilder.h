#pragma once
#include "Builder.h"

class NZoneBuilder : public Builder {
protected:
	virtual std::vector<double> construct_coef_matrix() override final;
	virtual double calc_dndx(double x) override final;
	virtual double get_Band_gap(double x) override;
public:
	NZoneBuilder(std::vector<double> &d_alfa);
	virtual void integrate_continuity_eq() override;
	virtual void write_to_file(const char* filename) override final;
	virtual double get_majority_carriers(double x) override final;
};