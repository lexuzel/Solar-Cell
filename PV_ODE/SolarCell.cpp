#include "pch.h"

void SolarCell::setBuilder(Builder &b) {
	builder = &b;
}
void SolarCell::execute() {
	builder->integrate_continuity_eq();
	if (delta_alfa == 0.0) {
	delta_alfa = builder->calc_delta();
	}

}
void SolarCell::writeToFile(const char* address) {
	builder->write_to_file(address);
}
double SolarCell::getDalfa() {
	return delta_alfa;
}