#pragma once
#include "Builder.h"

class SolarCell {
private:
	std::vector<double> delta_alfa;
	Builder* builder;
public:
	double ph_current = 0.0;
	void setBuilder(Builder &b);
	void execute();
	void writeToFile(const char* address);
	std::vector<double> getDalfa();
};
