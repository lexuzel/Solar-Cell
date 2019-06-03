// PV_ODE.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

int main(){
	SolarCell solar_cell;
	auto builder1 = GradedPZoneBuilder();
	solar_cell.setBuilder(builder1);
	solar_cell.execute();
	solar_cell.writeToFile("pzone.txt");

	auto delta_alfa = solar_cell.getDalfa();

	auto builder2 = PNJunctionBuilder(delta_alfa);
	solar_cell.setBuilder(builder2);
	solar_cell.execute();
	solar_cell.writeToFile("pn_junc.txt");

	delta_alfa = solar_cell.getDalfa();

	auto builder3 = NZoneBuilder(delta_alfa);
	solar_cell.setBuilder(builder3);
	solar_cell.execute();
	solar_cell.writeToFile("nzone.txt");

	std::cout << "Photo-current = " << solar_cell.ph_current;
}

