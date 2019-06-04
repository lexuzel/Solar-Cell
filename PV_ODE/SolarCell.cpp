#include "pch.h"

void SolarCell::setBuilder(Builder &b) {
	builder = &b;
}
void SolarCell::execute() {
	builder->calc_photo_carriers();
	ph_current += builder->photo_current;
}
void SolarCell::writeToFile(const char* address) {
	builder->write_to_file(address);
}
std::vector<double> SolarCell::getDalfa() {
	return delta_alfa;
}