#include "pch.h"

void SolarCell::setBuilder(Builder &b) {
	builder = &b;
}
void SolarCell::execute() {
	builder->integrate_continuity_eq();
	delta_alfa = builder->integrate_abs_coef(0, width_n);
}
void SolarCell::writeToFile(const char* address) {
	builder->write_to_file(address);
}
double SolarCell::getDalfa() {
	return delta_alfa;
}