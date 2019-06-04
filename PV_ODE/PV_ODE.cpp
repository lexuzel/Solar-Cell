// PV_ODE.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

int main(){
	SolarCell solar_cell;
	auto builder1 = EntireZoneBuilder();
	solar_cell.setBuilder(builder1);
	solar_cell.execute();
	solar_cell.writeToFile("entire-zone.txt");

	std::cout << "Photo-current = " << solar_cell.ph_current;
}

