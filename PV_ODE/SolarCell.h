#pragma once
#include "Builder.h"

class SolarCell {
private:
	Builder builder;
public:
	void setBuilder(const Builder &b);
	void execute();
	void writeToFile(const char* address);
};
