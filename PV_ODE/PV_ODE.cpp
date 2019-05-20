// PV_ODE.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

int main(){
	SolarCell solar_cell;
	auto builder1 = PZoneBuilder();
	solar_cell.setBuilder(builder1);
	solar_cell.execute();
	solar_cell.writeToFile("pzone.txt");

//	double delta_alfa = solar_cell.getDalfa();

	auto builder2 = GradedPZoneBuilder();
	solar_cell.setBuilder(builder2);
	solar_cell.execute();
	solar_cell.writeToFile("nzone.txt");
}

