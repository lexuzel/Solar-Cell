#pragma once
#include "Builder.h"
#include <map>

class SolarCell {
private:
	Builder* builder;
public:
	std::vector<std::array<double, 5>> summary;
	double ph_current = 0.0;
	double rec_current = 0.0;
	double minor_current = 0.0;
	double equ_current = 0.0;

	void setBuilder(Builder &b);
	void execute(double Eg_0 = 1.5);
	void execute_all(double begin, double end, double step);
	void writeToFile(const char* address);
	void writeSummary(const char* address);
};
