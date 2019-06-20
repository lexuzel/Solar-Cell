#pragma once
#include "Builder.h"

class EntireZoneBuilder : public Builder {
	std::vector<double> electrons, holes;
	std::vector<double> eq_electrons, eq_holes;
	int min_index;
protected:
	virtual double get_Eq(double x) override final;
	virtual std::pair<double, double> get_Energy_level(double x) override final;

	virtual double get_photoCurrent() override final;
	virtual double get_recombCurrent(std::vector<double> &table) override final;
	virtual double get_minorCurrent(std::vector<double> &table) override final;
	virtual double get_quantum_eff(std::vector<double>& table) override final;
public:
	EntireZoneBuilder(double width_p = 5.0e-6, double width_n = 2.0e-5);
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers(double voltage) override final;
};