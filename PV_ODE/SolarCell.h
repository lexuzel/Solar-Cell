#pragma once
#include "Builder.h"
#include <map>

class SolarCell {
private:
	Builder* builder;
public:
	std::vector<std::array<double, 8>> summary;
	double ph_current = 0.0;
	double minor_current = 0.0;
	double rec_current = 0.0;
	double equ_current = 0.0;

	void setBuilder(Builder &b);
	void execute(double width_p, double width_n, double Eg0 = 1.5, double S = 1.0e4);

	void investigateEg(double width_p, double width_n, double S, double begin, double end, double step);
	void investigateS(double width_p, double width_n, double Eg0, double begin, double end, double step);
	void investigateWidth(double begin, double end, double step, double width_n, double Eg0, double S);

	void writeToFile(const char* address);
	void writeSummary(const char* address);
};
