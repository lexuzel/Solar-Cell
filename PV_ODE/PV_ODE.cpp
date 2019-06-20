// PV_ODE.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

int main(){
	using clock = std::chrono::steady_clock;
	using ms = std::chrono::milliseconds;

	auto point1 = clock::now();

	SolarCell solar_cell;
	EntireZoneBuilder builder(4.5e-6);
	solar_cell.setBuilder(builder);
//	solar_cell.execute();

	solar_cell.execute_all(1.72, 1.82, 0.01);
	solar_cell.writeSummary("summary.txt");
//	solar_cell.writeToFile("entire-zone.txt");

	//double correction = 1.0e-3;

	//std::cout << " Photo-current = " << solar_cell.ph_current / 10 << " [mA / cm2]" << "\n";
	//std::cout << " Minority-current = " << solar_cell.minor_current / 10 << " [mA / cm2]" << "\n";
	//std::cout << " Shade current = " << solar_cell.equ_current / 10 * correction << " [mA / cm2]" << "\n";

	//std::ofstream fout;
	//fout.open("vax.txt");

	//double photo_cur = solar_cell.ph_current + solar_cell.minor_current;
	//double Voc = k_b * T / q_e * log(photo_cur / (solar_cell.equ_current * correction) + 1);

	//std::cout << " Short-curcuit current = " << photo_cur / 10 << " [mA / cm2]" << "\n";
	//std::cout << " Open-curcuit voltage = " << Voc << " [V]" << "\n";

	//double max_power{ 0 }, max_current, max_voltage;
	//fout.width(15);
	//fout << "Voltage \t Current \t Power \n\n";
	//for (double volt = 0.0; volt < Voc + 0.005; volt += 0.005) {
	//	double current = (photo_cur - solar_cell.equ_current * correction * (exp(q_e * volt / k_b / T) - 1)) / 10;
	//	double power = current * volt;
	//	fout.scientific;
	//	fout.width(15);
	//	fout << volt << "\t" << current << "\t" << power << "\n";
	//	if (power > max_power) {
	//		max_power = power;
	//		max_current = current;
	//		max_voltage = volt;
	//	}
	//}
	//fout.close();

	//std::cout << std::dec;
	//std::cout << "\n MPPT-power = " << max_power << " [mW / cm2]" << "\n";
	//std::cout << " MPPT-current = " << max_current << " [mA / cm2]" << "\n";
	//std::cout << " MPPT-voltage = " << max_voltage << " [V]" << "\n";
	//double FF = max_power / (photo_cur/10) / Voc;
	//std::cout << " Fill Factor = " << FF << "\n";
	//double efficiency = (photo_cur/10) * Voc * FF / 100;
	//std::cout << " Efficiency = " << efficiency * 100 << " [%]\n";

	std::cout << "\nExecution time: " << std::chrono::duration_cast<ms>(clock::now() - point1).count() << "\n";
}

