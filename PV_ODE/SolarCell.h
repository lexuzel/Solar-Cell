#pragma once
#include "Builder.h"

class SolarCell {
private:
	Builder* builder;
public:
	void setBuilder(Builder &b);
	void execute();
	void writeToFile(const char* address);
};
