#pragma once

/* ����, � ����� �������� ��������� �������� */

/* Բ���Ͳ ��������� */

const double pi = 3.14159265359;
const double plank_h = 6.626E-34;			// ����� ������		
const long speed_of_light = 299792458;		// �������� �����
const double q_e = 1.6021766208E-19;		// ����� ��������� (������)
const double m0 = 9.1093826E-31;			// ���� ������� ���������
const double k_b = 1.38064852E-23;			// ����� ���������
const double T = 300.0;						// ����������� � ��������

// ���������� �������� ����� [m]
const double x_1 = 3.5E-6;					// ʳ���� �������� �- ������ - ������� ��� (n-������� �������)
const double x_0 = x_1 + 1.5E-6;			// ������� p-n �������� ��������
const double x_2 = x_0 + 1.5E-6;			// ʳ���� ��� (�-������� �������) - ������� ��������� n- ������
const double x_d = x_2 + 1.0E-5;			// ʳ���� ��������� n- ������

// ������������ ������ [1/m3]
const double NA = 2.0E18;					// ������������ ��������� � �-������
const double ND = 2.0E18;					// ������������ ������ � n-������

// ���������� ���������� ��� Eg = hv
const double alfa0 = 1.0E7;					// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>

// ��������� ����糿 ��������� � ���� [m2/(V*s)]
const double D_diff_n = 1.0E-3;				// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>
const double D_diff_p = 1.0E-4;				// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>

// ��������� ��������� � ���� [m/s * m/V]
const double mu_n = 0.04;					// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>
const double mu_p = 0.004;					// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>

// ��� ����� ��������� � ���� 
const double tau_n = 1.0E-9;				// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>
const double tau_p = 1.0E-9;				// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>

// �������� ���������� ����������� [m/s]
const double S = 1.0E4;						// Optimization of graded band gap CdHgTe solar cells <A.Bouazzi>

// ������ ���������� ���� (������ �����) [eV]
const double Eg_0 = 2.76;
const double Eg_x1 = 2.0;
const double Eg_x2 = 1.5;

// �������� ���� ��������� � ���� [.]
const double me_p = 0.064 * Eg_x1;
const double mh_p = 0.35;
const double me_n = 0.1;
const double mh_n = 0.016;

// ������ ��������� �������� [eV]
//const double bias = (Eg_x2 - Eg_x1) / 2 + 0.75 * k_b * T * log(me_p * mh_n / mh_p / me_n);
const double bias = 1.0;

// г���� ���� [eV]
//const double fermi_level = Eg_x1 / 2 - 0.75 * k_b * T * log(me_p / mh_p);
const double fermi_level = 0.25;


/* ��������� ��������� ������ */

const long K_x = 5000;						// ѳ��� ������������� ��� � 

// ���� ���� ��� �

const double step_x_n = x_1 / K_x;	
const double step_x_p = (x_d - x_2) / K_x;
const double step_x_pn = (x_2 - x_1) / K_x;
