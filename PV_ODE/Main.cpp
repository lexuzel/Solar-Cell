
#include "pch.h"

int main(){
	Cell_Intrinsic n_zone(0.0, 0.0);
	n_zone.integrate_continuity_eq();
	n_zone.write_to_file("N_Table.txt");

//	double delta_alfa = n_zone.integrate_abs_coef(0, width_n);

//	Cell_Intrinsic p_zone(2.0, delta_alfa);
//	p_zone.integrate_continuity_eq();
//	p_zone.write_to_file("P_Table.txt");

	return 0;
}

