#include "pch.h"

void SolarCell::setBuilder(const Builder &b) {
	builder = b;
}
void SolarCell::execute() {
	builder.integrate_continuity_eq();
}
void SolarCell::writeToFile(const char* address) {
	builder.write_to_file(address);
}