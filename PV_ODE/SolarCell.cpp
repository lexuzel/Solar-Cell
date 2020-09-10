#include "pch.h"

void SolarCell::setBuilder(Builder &b) {
	builder = &b;
}

void SolarCell::execute(double width_p, double width_n, double Eg0, double S) {
	builder->calc_photo_carriers(width_p, width_n, Eg0, S);
	ph_current = builder->photo_current;
	minor_current = builder->minor_current;
	rec_current = builder->rec_current;
	equ_current = builder->equ_current;

	summary.push_back({ width_p, Eg0, builder->maxEg0, S, builder->photo_current, builder->minor_current, builder->rec_current, builder->equ_current });
}

void SolarCell::investigateEg(double width_p, double width_n, double S, double begin, double end, double step){
	for (double eg = begin; eg < end + step / 2; eg += step) {
		execute(width_p, width_n, eg, S);
	}
}

void SolarCell::investigateS(double width_p, double width_n, double Eg0, double begin, double end, double step) {
	for (double S = begin; S < end + step / 2; S *= step) {
		execute(width_p, width_n, Eg0, S);
	}
}

void SolarCell::investigateWidth(double begin, double end, double step, double width_n, double Eg0, double S)
{
	for (double width_p = begin; width_p < end + step / 2; width_p += step) {
		execute(width_p, width_n, Eg0, S);
	}
}

void SolarCell::writeToFile(const char* address) {
	builder->write_to_file(address);
}

enum Summary { WP, EG, EG0MAX, S, PH_CURRENT, MINOR_CURRENT, REC_CURRENT, EQU_CURRENT };

void SolarCell::writeSummary(const char * address){
	std::ofstream fout;
	fout.open(address);

	fout << "PZone width \t Eg_0 \t Eg0_max \t S \t photo current \t minor current \t Short-curcuit current \t Shade current \t Open-curcuit voltage \t MPPT power \t MPPT current \t MPPT voltage \t Fill Factor \t Efficiency \n\n";
	for (const auto &a : summary) {
		double photo_cur = a[MINOR_CURRENT] + a[REC_CURRENT];
		photo_cur = photo_cur < 0.0 ? 0.0 : photo_cur;
		double Voc = k_b * T / q_e * log(photo_cur / a[EQU_CURRENT] + 1);

		double max_power{ 0 }, max_current, max_voltage;
		for (double volt = 0.0; volt < Voc + 0.005; volt += 0.005) {
			double current = (photo_cur - a[EQU_CURRENT] * (exp(q_e * volt / k_b / T) - 1)) / 10;
			double power = current * volt;
			if (power > max_power) {
				max_power = power;
				max_current = current;
				max_voltage = volt;
			}
		}
		double FF = photo_cur != 0.0 ? max_power / (photo_cur) / Voc : 0.0;
		double efficiency = (photo_cur) * Voc * FF / 100;
		fout << a[WP] << "\t" 
			<< a[EG] << "\t"
			<< a[EG0MAX] << "\t"
			<< a[S] * 100 << "\t" 
			<< a[REC_CURRENT] / 10 << "\t"
			<< a[MINOR_CURRENT] / 10 << "\t"
			<< photo_cur / 10 << "\t" 
			<< a[EQU_CURRENT] / 10 << "\t" 
			<< Voc << "\t"
			<< max_power << "\t" 
			<< max_current << "\t" 
			<< max_voltage << "\t" 
			<< FF << "\t" 
			<< efficiency * 100<< "\n";
	}
	fout.close();
}
