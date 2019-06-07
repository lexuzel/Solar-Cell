#include "pch.h"

void SolarCell::setBuilder(Builder &b) {
	builder = &b;
}

void SolarCell::execute(double Eg_0) {
	builder->calc_photo_carriers(Eg_0);
	ph_current = builder->photo_current;
	rec_current = builder->recomb_current;
	minor_current = builder->minor_current;
	equ_current = builder->equ_current;

	summary.push_back({ Eg_0, ph_current/10, rec_current/10, minor_current/10, equ_current/10 });

	builder->photo_current = 0.0;
	builder->recomb_current = 0.0;
	builder->minor_current = 0.0;
	builder->equ_current = 0.0;
}

void SolarCell::execute_all(double begin, double end, double step){
	for (double Eg_0 = begin; Eg_0 < end; Eg_0 += step) {
		execute(Eg_0);
	}
}

void SolarCell::writeToFile(const char* address) {
	builder->write_to_file(address);
}

void SolarCell::writeSummary(const char * address){
	std::ofstream fout;
	fout.open(address);

	fout << "Eg_0 \t Short-curcuit current \t Open-curcuit voltage \t MPPT power \t MPPT current \t MPPT voltage \t Fill Factor \t Efficiency \n\n";
	for (const auto &a : summary) {
		double photo_cur = a[1] + a[3];
		double Voc = k_b * T / q_e * log(photo_cur / a[4] + 1);

		double max_power{ 0 }, max_current, max_voltage;
		for (double volt = 0.0; volt < Voc + 0.005; volt += 0.005) {
			double current = (photo_cur - a[4] * (exp(q_e * volt / k_b / T) - 1));
			double power = current * volt;
			if (power > max_power) {
				max_power = power;
				max_current = current;
				max_voltage = volt;
			}
		}
		double FF = max_power / (photo_cur / 10) / Voc;
		double efficiency = (photo_cur / 10) * Voc * FF / 100;
		fout << a[0] << "\t" << photo_cur << "\t" << Voc << "\t" 
			<< max_power << "\t" << max_current << "\t" << max_voltage << "\t" 
			<< FF << "\t" << efficiency * 100<< "\n";
	}
	fout.close();
}
