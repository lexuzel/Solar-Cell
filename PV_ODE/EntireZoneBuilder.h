#pragma once
#include "Builder.h"

class EntireZoneBuilder : public Builder {
	std::vector<double> electrons, holes;
	std::vector<double> eq_electrons, eq_holes;
	
protected:
	virtual double get_Eq(double x) override final;
	virtual std::pair<double, double> get_Energy_level(double x) override final;
	void getDepletionRegionCoord();

	enum Zone { PZONE, IZONE, NZONE };
	double intrCarrierConc(Zone zone);

	virtual double get_photoCurrent() override final;
	virtual double get_minorCurrent() override final;
	virtual double get_recCurrent(std::vector<double>& table) override final;

	double getAverageField(double p1, double p2);
	double get_satCurrent();

public:
	EntireZoneBuilder();
	virtual void write_to_file(const char* filename) override final;
	virtual void calc_photo_carriers(double width_p, double width_n, double Eg0, double S) override final;
};