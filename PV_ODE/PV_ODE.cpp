// PV_ODE.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"

int main(){
	SolarCell solar_cell;
	Builder builder = GradedPZoneBuilder();
	solar_cell.setBuilder(builder);
	solar_cell.execute();
	solar_cell.writeToFile("test.txt");
}

