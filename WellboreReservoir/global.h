#pragma once
#include<vector>

using namespace std;


#define STBtoSCF 5.615
#define Psc 101325.01
#define Tsc 288.15
#define tol1 1e-8
#define Rgas 8.3144621
#define piNum 3.141592653589793
#define Tref_ 298.15
#define g_acc 9.81	//	m / s

static bool AlmostEqual(double Expected, double Actual, double tol) {
	if (abs(Expected - Actual) < tol)
	{
		return true;
	}
	else {
		return false;
	}
}

static double findMin(double* x, int N) {
	register int i;
	double min = x[0];
	for (i = 1; i < N; i++) {
		if (min > x[i]) {
			min = x[i];
		}
	}
	return min;
}

static double findMax(double* x, int N) {
	register int i;
	double max = x[0];
	for (i = 1; i < N; i++) {
		if (max < x[i]) {
			max = x[i];
		}
	}
	return max;

}

static void normalizArr(double* x, int N) {
	double sum = 0.0;
	register int i;
	for (i = 0; i < N; i++) {
		sum += x[i];
	}
	for (i = 0; i < N; i++) {
		x[i] /= sum;
	}
}


enum fluidPhase{Total = 0, Liquid = 1, Gas = 2};

struct pvt_dataStruct
{
	double P;
	double T;
	int Np;
	vector<double> comp;
	vector<double> phi;

};


static double F_to_C(double F) {
	return (F - 32.0)*5.0 / 9.0;
}
static double F_to_K(double F) {
	return (F - 32.0)*5.0 / 9.0 + 273.15;
}
static double F_to_R(double F) {
	return (F + 459.67);
}

static double C_to_F(double C) {
	return C * 9.0 / 5.0 + 32.0;
}

static double C_to_K(double C) {
	return C + 273.15;
}

static double C_to_R(double C) {
	return C * 9.0 / 5.0 + 491.67;
}

static double R_to_C(double R) {
	return (R - 491.67) * 5.0 / 9.0;
}

static double R_to_K(double R) {
	return R * 5.0 / 9.0;
}

static double R_to_F(double F) {
	return (F - 459.67);
}

static double K_to_C(double K){
	return (K - 273.15);
}

static double K_to_F(double K) {
	return (K - 273.15)*9.0 / 5.0 + 32.0;
}

static double K_to_R(double K) {
	return (K * 1.8);
}

static double m_to_foot(double m) {
	return (m * 3.280839895);
}

static double m_to_inch(double m) {
	return (m * 39.37007874);
}

static double inch_to_meter(double inch) {
	return (inch * 0.0254);
}

static double inch_to_foot(double inch) {
	return (inch * 0.0833333333);
}

static double foot_to_meter(double foot) {
	return (foot * 0.3048);
}

static double foot_to_inch(double foot) {
	return (foot * 12.0);
}

static double bar_to_atm(double bar) {
	return (bar *  0.98692316931);
}

static double bar_to_Pas(double bar) {
	return (bar * 1e+5);
}

static double bar_to_psi(double bar) {
	return (bar * 14.503773801);
}

static double atm_to_bar(double atm) {
	return (atm * 1.0132501);
}

static double atm_to_Pas(double atm) {
	return (atm * 101325.01);
}

static double atm_to_psi(double atm) {
	return (atm * 14.695950254);
}

static double Pas_to_atm(double Pas) {
	return (Pas * 0.0000098692316931);
}

static double Pas_to_bar(double Pas) {
	return (Pas * 1e-5);
}

static double Pas_to_psi(double Pas) {
	return (Pas * 0.00014503773801);
}

static double psi_to_atm(double psi) {
	return (psi * 0.068045957064);
}

static double psi_to_bar(double psi) {
	return (psi * 0.0689475728);
}

static double psi_to_Pas(double psi) {
	return (psi * 6894.75728);
}

static double Day_to_s(double Day) {
	return (Day * 86400.0);
}

static double s_to_Day(double s) {
	return (s * 1.1574074 * 1e-5);
}

static double scf_to_scm(double scf) {
	return (scf * 0.028316846592);
}

static double scf_to_STB(double scf) {
	return (scf * 0.178107606679035);
}

static double scm_to_scf(double scm) {
	return scm * 35.3146667214886;
}

static double scm_to_gallon(double scm) {
	return (scm * 264.172052358148);
}

static double STB_to_scf(double STB) {
	return (STB * 5.61458333333333);
}

static double STB_to_scm(double STB) {
	return (STB * 0.158987294928);
}

static double STB_to_gallon(double STB) {
	return (STB * 42.0);
}

static double scfPerSTB_to_scmPerscm(double scfPerSTB) {
	return (scfPerSTB * 0.178107606679035);
}

static double scmPerscm_To_scfPerSTB(double scmPerscm) {
	return (scmPerscm * 5.61458333333333);
}

static double PaS_to_cp(double Pas) {
	return Pas * 1000.0;
}

static double cp_to_PaS(double cp) {
	return cp / 1000.0;
}

static double bmft3_to_kgm3(double bmft3) {
	return bmft3 / 0.062428;
}

static double kgm3_to_bmft3(double kgm3) {
	return kgm3 / 0.062428;
}

static double dynePerCen_to_NPerM(double dynePerCen) {
	return 0.001 * dynePerCen;
}

static double NPerM_to_dynePerCen(double NPerM) {
	return NPerM * 1000.0;
}