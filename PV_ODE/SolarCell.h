#pragma once
#include "Builder.h"

class SolarCell {
private:
	double delta_alfa;
	Builder* builder;
public:
	void setBuilder(Builder &b);
	void execute();
	void writeToFile(const char* address);
	double getDalfa();
};
