#include "BlackOil.h"



BlackOil::BlackOil()
{

}


BlackOil::~BlackOil()
{
}

BlackOil::BlackOil(double T, double P, double gammaGas, double API, double Rs, double Co,
	string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
	string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod) {
	if ((units == "Metric") || (units == "metric") || (units == "SI")) {
		this->T = T;
		this->P = P;
		this->Rs = Rs;
		this->Co = Co;
	}
	else if ((units == "Field") || (units == "field")) {
		this->T = F_to_K(T);
		this->P = psi_to_Pas(P);
		this->Rs = scfPerSTB_to_scmPerscm(Rs);
		this->Co = Co / psi_to_Pas(1.0);
	}
	this->gammaGas = gammaGas;
	this->API = API;
	gammaOil = gammaOilCal();
	Pseudocritical_props(Pseudocritical_method);
	MwOil = stockTankMw(MwMethod);
	MwGas = 28.97 * gammaGas / 1000.0;
	this->muGasMethod = muGasMethod;
	this->muDeadOilMethod = muDeadOilMethod;
	this->muSatOilMethod = muSatOilMethod;
	this->ZFactorMethod = ZFactorMethod;
	this->BobMethod = BobMethod;
}


BlackOil::BlackOil(double gammaGas, double API, double Rs, double Co,
	string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
	string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod) {

	this->gammaGas = gammaGas;
	this->API = API;
	gammaOil = gammaOilCal();
	Pseudocritical_props(Pseudocritical_method);
	MwOil = stockTankMw(MwMethod);
	MwGas = 28.97 / 1000.0 * gammaGas;
	this->muGasMethod = muGasMethod;
	this->muDeadOilMethod = muDeadOilMethod;
	this->muSatOilMethod = muSatOilMethod;
	this->satOilVisMethod = satOilVisMethod;
	this->ZFactorMethod = ZFactorMethod;
	this->BobMethod = BobMethod;
	if ((units == "Metric") || (units == "metric") || (units == "SI")) {
		this->Rs = Rs;
		this->Co = Co;
	}
	else if ((units == "Field") || (units == "field")) {
		this->Rs = scfPerSTB_to_scmPerscm(Rs);
		this->Co = Co / psi_to_Pas(1.0);
	}
}


BlackOil::BlackOil(double gammaGas, double API, double Rs, double Co, double waterSal,
	string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
	string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod) {

	this->waterSal = waterSal;
	this->gammaGas = gammaGas;
	this->API = API;
	gammaOil = gammaOilCal();
	Pseudocritical_props(Pseudocritical_method);
	MwOil = stockTankMw(MwMethod);
	MwGas = 28.97 / 1000.0 * gammaGas;
	this->muGasMethod = muGasMethod;
	this->muDeadOilMethod = muDeadOilMethod;
	this->muSatOilMethod = muSatOilMethod;
	this->satOilVisMethod = satOilVisMethod;
	this->ZFactorMethod = ZFactorMethod;
	this->BobMethod = BobMethod;
	if ((units == "Metric") || (units == "metric") || (units == "SI")) {
		this->Rs = Rs;
		this->Co = Co;
	}
	else if ((units == "Field") || (units == "field")) {
		this->Rs = scfPerSTB_to_scmPerscm(Rs);
		this->Co = Co / psi_to_Pas(1.0);
	}
}








double BlackOil::bubblePointPressure(string corr) {

	if (corr == "Standing") {
		/*The Standing (1947) correlation is based on a laboratory
			study of 22 different crude oils from California.
			It is best suited for oils with gravity less than 15 API.*/
		Pb = PbStanding();
	}
	else if (corr == "Lasater") {
		Pb = PbLasater();
	}
	else if (corr == "Glaso") {
		Pb = PbGlaso();
	}
	else if (corr == "Beggs-Vazquez") {
		Pb = PbBeggs_Vazquez();
	}
	else if (corr == "Dindoruk-Christman") {
		Pb = PbDindoruk_Christman();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate bubble point pressure.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return Pb;

}

double BlackOil::SaturatedGasOilRatio(string corr) {

	if (corr == "Standing") {
		Rs = RsStanding();
	}
	else if (corr == "Lasater") {
		Rs = RsLasater();
	}
	else if (corr == "Vazquez-Beggs") {
		Rs = RsVazquez_Beggs();
	}
	else if (corr == "Glaso") {
		Rs = RsGlaso();
	}
	else if (corr == "Dindoruk-Christman") {
		Rs = RsDindoruk_Christman();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate saturated Gas/Oil ratio.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return Rs;
}

double BlackOil::OFVF(string corr) {

	if (corr == "Standing") {
		Bob = BobStanding();
	}
	else if (corr == "Glaso") {
		Bob = BobGlaso();
	}
	else if (corr == "AlMarhoun") {
		Bob = BobAlMarhoun();
	}
	else if (corr == "Vazquez-Beggs") {
		Bob = BobVazquez_Beggs();
	}
	else if (corr == "Dindoruk-Christman") {
		Bob = BobDindoruk_Christman();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate oil formation volume factor.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return Bob;
}

void BlackOil::DeadOilViscosity(string corr) {
	if (corr == "Beal-Standing") {
		muDeadBeal_Standing();
	}
	else if (corr == "Beggs-Robinson") {
		muDeadBeggs_Robinson();
	}
	else if (corr == "Glaso") {
		muDeadGlaso();
	}
	else if (corr == "AlKhafaji") {
		muDeadAlKhafaji();
	}
	else if (corr == "Dindoruk-Christman") {
		muDeadDindoruk_Christman();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate dead oil viscosity.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}

}

double BlackOil::SatOilViscosity() {

	if (muSatOilMethod == "Standing") {
		return muSatStanding();
	}
	else if (muSatOilMethod == "Beggs-Robinson") {
		return muSatBeggs_Robinson();
	}
	else if (muSatOilMethod == "Bergman") {
		return muSatBergman();
	}
	else if (muSatOilMethod == "AlKhafaji") {
		return muSatAlKhafaji();
	}
	else if (muSatOilMethod == "Aziz") {
		return muSatAziz();
	}
	else if (muSatOilMethod == "Dindoruk-Christman") {
		return muSatDindoruk_Christman();
	}
	else {
		string msg;
		msg = "Undefined function: " + satOilVisMethod + " to caculate saturated oil viscosity.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
}

double BlackOil::GFVF() {
	Bg = BgCalculation();
	return Bg;
}

double BlackOil::GasViscosity(string corr) {
	if (corr == "Dempsey") {
		muGas = muGasDempsey();
	}
	else if (corr == "Lee-Gunzalez") {
		muGas = muGasLee_Gunzalez();
	}
	else if (corr == "Lucas") {
		muGas = muGasLucas();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate gas viscosity.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return muGas;
}

double BlackOil::ZFactor() {
	if (ZFactorMethod == "PVTsim") {
		Z = ZFactor_PVTsim();
	}
	else {
		string msg;
		msg = "Undefined function: " + ZFactorMethod + " to caculate gas Z factor (compressibility factor).";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return Z;
}

void BlackOil::Pseudocritical_props(string corr) {
	if (corr == "Sutton") {
		Pc_Tc_Sutton();
	}
	else if (corr == "T2") {
		gasPpc = psi_to_Pas(756.8 - 131.0*gammaGas - 3.6*gammaGas*gammaGas);
		gasTpc = R_to_K(169.2 + 349.5*gammaGas - 74.0*gammaGas*gammaGas);
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate pseudocritical properties.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
}

double BlackOil::PbStanding() {


	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempT = K_to_F(T);
	double A = pow((tempRs / gammaGas), 0.83) * pow(10.0, 0.00091*tempT - 0.0125 * API);
	Pb = psi_to_Pas(18.2 * (A - 1.4));
	return Pb;
}

double BlackOil::PbLasater() {

	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempT = K_to_R(T);

	double yg = (tempRs / 379.3) / (tempRs / 379.3 + 350 * gammaOil / MwOil);
	double A;

	if (gammaGas <= 0.6) {
		A = (0.679 * exp(2.786 * yg) - 0.323);
	}
	else {
		A = 8.26 * pow(yg, 3.56) + 1.95;
	}

	Pb = (tempT / gammaGas) * A;
	Pb = psi_to_Pas(Pb);
	return Pb;
}

double BlackOil::PbGlaso() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempT = K_to_F(T);
	double A = pow(tempRs / gammaGas, 0.816) * (pow(tempT, 0.172) / pow(API, 0.989));
	double temp = 1.7669 + 1.7447 * log10(A) - 0.30218*log10(A)*log10(A);
	Pb = pow(10, temp);

	Pb = psi_to_Pas(Pb);
	return Pb;

}

double BlackOil::PbBeggs_Vazquez() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempT = K_to_F(T);

	if (API <= 30.0) {
		double C1 = 0.0362; double C2 = 1.0937; double C3 = 25.724;
		double temp = tempRs / (C1 * gammaGas * exp(C3 * (API / (tempT + 459.67))));
		Pb = pow(temp, 1 / C2);
		// PVTsim formula
		//double temp = 27.64 * (tempRs / gammaGas) * pow(10.0, (-11.172*API / (tempT + 459.67)));
		//Pb = pow(temp, 0.9143);
	}
	else {
		double C1 = 0.0178; double C2 = 1.1870; double C3 = 23.9310;
		double temp = tempRs / (C1 * gammaGas * exp(C3 * (API / (tempT + 459.67))));
		Pb = pow(temp, 1 / C2);
		// PVTsim formula
		//double temp = 56.06 * (tempRs / gammaGas) * pow(10.0, (-10.393*API) / (tempT + 459.67));
		//Pb = pow(temp, 0.8425);
	}

	Pb = psi_to_Pas(Pb);
	return Pb;
}

double BlackOil::PbDindoruk_Christman() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempT = K_to_F(T);
	double a1 = 1.42828 * 1e-10; double a2 = 2.844591797; double a3 = -6.74896 * 1e-4;
	double a4 = 1.225226436; double a5 = 0.033383304; double a6 = -0.272945957; double a7 = -0.084226069;
	double a8 = 1.869979257; double a9 = 1.221486524; double a10 = 1.370508349; double a11 = 0.011688308;
	double temp1 = a1 * pow(tempT, a2) + a3 * pow(API, a4);
	double temp2 = (a5 + 2 * pow(tempRs, a6) / pow(gammaGas, a7)) * (a5 + 2 * pow(tempRs, a6) / pow(gammaGas, a7));
	double A = temp1 / temp2;
	Pb = a8 * ((pow(tempRs, a9) / pow(gammaGas, a10)) * pow(10, A) + a11);
	Pb = psi_to_Pas(Pb);
	return Pb;

}

double BlackOil::RsStanding() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double temp = (tempP / 18.2 + 1.4) * pow(10.0, 0.0125 * API - 0.00091 * tempT);
	//double temp = ((1.0/18.0) * tempP) * pow(10.0, 0.0125 * API) / pow(10, 0.00091*tempT);

	Rs = gammaGas * pow(temp, 1.205);
	Rs = scfPerSTB_to_scmPerscm(Rs);

	return Rs;
}

double BlackOil::RsLasater() {
	double tempT = K_to_R(T);
	double tempP = Pas_to_psi(P);
	double yg;
	if (tempP * gammaGas / tempT < 3.29) {
		yg = 0.359 * log(1.473 * tempP * gammaGas / tempT + 0.476);
	}
	else {
		yg = pow((0.121*tempP*gammaGas / tempT - 0.236), 0.281);
	}
	Rs = 132755 * gammaOil * yg / (MwOil*(1 - yg));
	Rs = scfPerSTB_to_scmPerscm(Rs);

	return Rs;
}

double BlackOil::RsVazquez_Beggs() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double C1, C2, C3;

	if (API <= 30.0) {
		C1 = 0.0362; C2 = 1.0937; C3 = 25.724;
	}
	else {
		C1 = 0.0178; C2 = 1.1870; C3 = 23.931;
	}
	Rs = C1 * gammaGas * pow(tempP, C2) * exp((C3 * API) / (tempT + 495.67));
	Rs = scfPerSTB_to_scmPerscm(Rs);

	return Rs;
}

double BlackOil::RsGlaso() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double temp1 = (-1.7447 + sqrt(3.044 + 1.20872*(1.7669 - log10(tempP)))) / -0.60436;
	temp1 = pow(10, temp1);
	double temp2 = pow(tempT, 0.172) / pow(API, 0.989);
	Rs = gammaGas * pow(temp1 / temp2, 1.0 / 0.816);
	Rs = scfPerSTB_to_scmPerscm(Rs);

	return Rs;
}

double BlackOil::RsDindoruk_Christman() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double a1 = 4.86996 * 1e-6; double a2 = 5.730982539; double a3 = 9.92510 * 1e-3;
	double a4 = 1.776179364; double a5 = 44.25002680; double a6 = 2.702889206; double a7 = 0.744335673;
	double a8 = 3.359754970; double a9 = 28.10133245; double a10 = 1.579050160; double a11 = 0.928131344;
	double temp1 = a1 * pow(API, a2) + a3 * pow(tempT, a4);
	double temp2 = (a5 + 2 * pow(API, a6) / pow(tempP, a7))*(a5 + 2 * pow(API, a6) / pow(tempP, a7));
	double A = temp1 / temp2;
	Rs = pow(((tempP / a8 + a9)*pow(gammaGas, a10)*pow(10, A)), a11);

	Rs = scfPerSTB_to_scmPerscm(Rs);

	return Rs;
}

double BlackOil::BobStanding() {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);

	double A = tempRs * sqrt(gammaGas / gammaOil) + 1.25*tempT;
	//Bob = 0.9759 + 12 * 1e-5 * pow(A, 1.2);
	Bob = 0.972 + 1.47 * 1e-4 * pow(A, 1.175);
	return Bob;
}

double BlackOil::BobStanding(double T, double Rs, double gammaGas, double gammaOil) {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);

	double A = tempRs * sqrt(gammaGas / gammaOil) + 1.25*tempT;
	Bob = 0.9759 + 12 * 1e-5 * pow(A, 1.2);

	return Bob;
}




double BlackOil::BobGlaso() {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double A = tempRs * pow(gammaGas / gammaOil, 0.526) + 0.968*tempT;
	double temp = -6.585 + 2.9133*log10(A) - 0.2768 * log10(A) * log10(A);
	Bob = pow(10.0, temp) + 1;

	return Bob;
}

double BlackOil::BobGlaso(double T, double Rs, double gammaGas, double gammaOil) {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double A = tempRs * pow(gammaGas / gammaOil, 0.526) + 0.968*tempT;
	double temp = -6.585 + 2.9133*log10(A) - 0.2768 * log10(A) * log10(A);
	Bob = pow(10.0, temp) + 1;

	return Bob;
}



double BlackOil::BobAlMarhoun() {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	Bob = 1.0 + (0.177342*1e-3)*tempRs + (0.220163*1e-3)*tempRs*gammaGas / gammaOil +
		(4.292580*1e-6)*tempRs*(tempT - 60.0)*(1.0 - gammaOil) + (0.528707*1e-3)*(tempT - 60.0);

	return Bob;
}

double BlackOil::BobAlMarhoun(double T, double Rs, double gammaGas, double gammaOil) {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	Bob = 1.0 + (0.177342*1e-3)*tempRs + (0.220163*1e-3)*tempRs*gammaGas / gammaOil +
		(4.292580*1e-6)*tempRs*(tempT - 60.0)*(1.0 - gammaOil) + (0.528707*1e-3)*(tempT - 60.0);

	return Bob;
}



double BlackOil::BobVazquez_Beggs() {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double C1, C2, C3;
	if (API <= 30.0) {
		C1 = 4.677 * 1e-4; C2 = 1.751*1e-5; C3 = -1.811*1e-8;
	}
	else {
		C1 = 4.670 * 1e-4; C2 = 1.100*1e-5; C3 = -1.337*1e-9;
	}
	Bob = 1.0 + C1 * tempRs + C2 * (tempT - 60.0)*(API / gammaGas) +
		C3 * tempRs * (tempT - 60.0) * (API / gammaGas);

	return Bob;
}

double BlackOil::BobVazquez_Beggs(double T, double Rs, double API, double gammaGas) {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double C1, C2, C3;
	if (API <= 30.0) {
		C1 = 4.677 * 1e-4; C2 = 1.751*1e-5; C3 = -1.811*1e-8;
	}
	else {
		C1 = 4.670 * 1e-4; C2 = 1.100*1e-5; C3 = -1.337*1e-9;
	}
	Bob = 1.0 + C1 * tempRs + C2 * (tempT - 60.0)*(API / gammaGas) +
		C3 * tempRs * (tempT - 60.0) * (API / gammaGas);

	return Bob;
}




double BlackOil::BobDindoruk_Christman() {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double a1 = 2.510755; double a2 = -4.852538; double a3 = 1.183500*1e+1; double a4 = 1.365428 * 1e+5;
	double a5 = 2.252880; double a6 = 1.007190 * 1e+1; double a7 = 4.450849 * 1e-1; double a8 = 5.352624;
	double a9 = -6.309052 * 1e-1; double a10 = 9.000749*1e-1; double a11 = 9.871766 * 1e-1;
	double a12 = 7.865146*1e-4; double a13 = 2.689173 * 1e-6; double a14 = 1.100001 * 1e-5;
	double temp1 = pow(tempRs, a1)*pow(gammaGas, a2) / pow(gammaOil, a3) + a4 * pow(tempT - 60.0, a5) + a6 * tempRs;
	double temp2 = (a8 + 2 * pow(tempRs, a9) / pow(gammaGas, a10) * (tempT - 60.0))
		*(a8 + 2 * pow(tempRs, a9) / pow(gammaGas, a10) * (tempT - 60.0));
	double A = pow(temp1, a7) / temp2;
	Bob = a11 + a12 * A + a13 * A * A + a14 * (tempT - 60.0) * API / gammaGas;
	return Bob;
}

double BlackOil::BobDindoruk_Christman(double T, double Rs, double gammaGas, double gammaOil, double API) {
	double tempT = K_to_F(T);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double a1 = 2.510755; double a2 = -4.852538; double a3 = 1.183500*1e+1; double a4 = 1.365428 * 1e+5;
	double a5 = 2.252880; double a6 = 1.007190 * 1e+1; double a7 = 4.450849 * 1e-1; double a8 = 5.352624;
	double a9 = -6.309052 * 1e-1; double a10 = 9.000749*1e-1; double a11 = 9.871766 * 1e-1;
	double a12 = 7.865146*1e-4; double a13 = 2.689173 * 1e-6; double a14 = 1.100001 * 1e-5;
	double temp1 = pow(tempRs, a1)*pow(gammaGas, a2) / pow(gammaOil, a3) + a4 * pow(tempT - 60.0, a5) + a6 * tempRs;
	double temp2 = (a8 + 2 * pow(tempRs, a9) / pow(gammaGas, a10) * (tempT - 60.0))
		*(a8 + 2 * pow(tempRs, a9) / pow(gammaGas, a10) * (tempT - 60.0));
	double A = pow(temp1, a7) / temp2;
	Bob = a11 + a12 * A + a13 * A * A + a14 * (tempT - 60.0) * API / gammaGas;
	return Bob;
}


double BlackOil::BobFunc() {
	if (BobMethod == "Standing") {
		Bob = BobStanding();
	}
	else if (BobMethod == "Glaso") {
		Bob = BobGlaso();
	}
	else if (BobMethod == "AlMarhoun") {
		Bob = BobAlMarhoun();
	}
	else if (BobMethod == "Vazquez-Beggs") {
		Bob = BobVazquez_Beggs();
	}
	else if (BobMethod == "Dindoruk-Christman") {
		Bob = BobDindoruk_Christman();
	}
	else {
		cout << "Wrong Method!" << endl;
	}
	return Bob;
}

double BlackOil::BoUnderSat() {
	BobFunc();
	BoUnd = Bob * exp(Co * (Pb - P));
	return BoUnd;
}

double BlackOil::BwFunc() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double dvWT = -1.0001 * 1e-2 + 1.33391 * 1e-4 * tempT + 5.50654 * 1e-7 * tempT * tempT;
	double dvWP = -1.95301 * 1e-9 * tempP * tempT - 1.72834 * 1e-13 * tempP * tempP * tempT -
		3.58922 * 1e-7 * tempP - 2.25341 * 1e-10 * tempP * tempP;

	Bw = (1.0 + dvWT) * (1.0 - dvWP);
	return Bw;
}

double BlackOil::waterDensityFunc() {
	rhoW = 62.368 + 0.438603 * waterSal + 1.60074 * 1e-3 * waterSal * waterSal;
	rhoW = bmft3_to_kgm3(rhoW);
	return rhoW;
}



double BlackOil::muDeadBeal_Standing() {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	muDeadOil = (0.32 + (1.8 * 1e+7) / pow(API, 4.53))*pow(360.0 / (tempT + 200.0), A);
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadBeal_Standing(double T, double API) {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	muDeadOil = (0.32 + (1.8 * 1e+7) / pow(API, 4.53))*pow(360.0 / (tempT + 200.0), A);
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}


double BlackOil::muDeadBeggs_Robinson() {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	if (tempT >= 70) {
		muDeadOil = -1 + pow(10, pow(tempT, -1.163) * exp(6.9824 - 0.04658 * API));
	}
	else {
		double A0 = 22.33 - 0.194*API + 0.00033 * API * API;
		double A1 = -3.20 + 0.0185 * API;
		muDeadOil = exp(A0 + A1 * log(tempT + 310.0)) - 1.0;
	}
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadBeggs_Robinson(double T, double API) {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	if (tempT >= 70) {
		muDeadOil = -1 + pow(10, pow(tempT, -1.163) * exp(6.9824 - 0.04658 * API));
	}
	else {
		double A0 = 22.33 - 0.194*API + 0.00033 * API * API;
		double A1 = -3.20 + 0.0185 * API;
		muDeadOil = exp(A0 + A1 * log(tempT + 310.0)) - 1.0;
	}
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}




double BlackOil::muDeadGlaso() {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	muDeadOil = (3.141*1e+10)*pow(tempT, -3.444)*pow(log10(API), (10.313*log10(tempT) - 36.447));
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadGlaso(double T, double API) {
	double tempT = K_to_F(T);
	double A = pow(10.0, 0.43 + (8.33 / API));
	muDeadOil = (3.141*1e+10)*pow(tempT, -3.444)*pow(log10(API), (10.313*log10(tempT) - 36.447));
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}



double BlackOil::muDeadAlKhafaji() {
	double tempT = K_to_F(T);
	muDeadOil = pow(10.0, 4.9563 - 0.00488*tempT) / pow(API + tempT / 30.0 - 14.29, 2.709);
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadAlKhafaji(double T, double API) {
	double tempT = K_to_F(T);
	muDeadOil = pow(10.0, 4.9563 - 0.00488*tempT) / pow(API + tempT / 30.0 - 14.29, 2.709);
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}


double BlackOil::muDeadDindoruk_Christman() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);

	double a1 = 14.505357625; double a2 = -44.868655416; double a3 = 9.36579 * 1e+9; double a4 = -4.194017808;
	double a5 = 3.1461171 * 1e-9; double a6 = 1.517652716; double a7 = 0.010433654; double a8 = -0.000776880;

	double A = a1 * log10(tempT) + a2;
	muDeadOil = a3 * pow(tempT, a4) * pow(log10(API), A) / (a5 * pow(tempP, a6) + a7 * pow(tempRs, a8));
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadDindoruk_Christman(double T, double P, double API) {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double tempRs = scmPerscm_To_scfPerSTB(Rs);

	double a1 = 14.505357625; double a2 = -44.868655416; double a3 = 9.36579 * 1e+9; double a4 = -4.194017808;
	double a5 = 3.1461171 * 1e-9; double a6 = 1.517652716; double a7 = 0.010433654; double a8 = -0.000776880;

	double A = a1 * log10(tempT) + a2;
	muDeadOil = a3 * pow(tempT, a4) * pow(log10(API), A) / (a5 * pow(tempP, a6) + a7 * pow(tempRs, a8));
	muDeadOil = cp_to_PaS(muDeadOil);
	return muDeadOil;
}

double BlackOil::muDeadOilFunc() {
	if (muDeadOilMethod == "Beal-Standing") {
		return muDeadBeal_Standing();
	}
	else if (muDeadOilMethod == "Beggs-Robinson") {
		return muDeadBeggs_Robinson();
	}
	else if (muDeadOilMethod == "Glaso") {
		return muDeadGlaso();
	}
	else if (muDeadOilMethod == "AlKhafaji") {
		return muDeadAlKhafaji();
	}
	else if (muDeadOilMethod == "Dindoruk-Christman") {
		return muDeadDindoruk_Christman();
	}
	else {
		cout << "Dead oil viscosity mrthod is wrong, Black Oil PVT******" << endl;
	}


}


double BlackOil::muSatStanding() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double temp1 = (-7.4*1e-4)*tempRs + (2.2*1e-7)*tempRs*tempRs;
	double A1 = pow(10.0, temp1);
	double A2 = 0.68 / pow(10.0, (8.62*1e-5)*tempRs) + 0.25 / pow(10.0, (1.1*1e-3)*tempRs) + 0.062 / pow(10.0, (3.74*1e-3)*tempRs);

	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}

double BlackOil::muSatStanding(double Rs, double muDeadOil) {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double temp1 = (-7.4*1e-4)*tempRs + (2.2*1e-7)*tempRs*tempRs;
	double A1 = pow(10.0, temp1);
	double A2 = 0.68 / pow(10.0, (8.62*1e-5)*tempRs) + 0.25 / pow(10.0, (1.1*1e-3)*tempRs) + 0.062 / pow(10.0, (3.74*1e-3)*tempRs);

	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}



double BlackOil::muSatBeggs_Robinson() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double A1 = 10.715*pow((tempRs + 100.0), -0.515);
	double A2 = 5.44 * pow((tempRs + 150.0), -0.338);
	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}

double BlackOil::muSatBergman() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double A1 = exp(4.768 - 0.8359*log(tempRs + 300.0));
	double A2 = 0.555 + 133.5 / (tempRs + 300.0);
	muSatOil = A1 * pow(muDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}


double BlackOil::muSatAziz() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double A1 = 0.20 + (0.80*pow(10.0, -0.00081*tempRs));
	double A2 = 0.43 + (0.57*pow(10.0, -0.00072*tempRs));
	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}

double BlackOil::muSatAlKhafaji() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double A0 = log10(tempRs);
	double A1 = 0.247 + 0.2824 * A0 + 0.5657 * A0 * A0 - 0.4065 * A0 * A0 *A0 + 0.0631* A0 * A0 * A0 * A0;
	double A2 = 0.894 + 0.0546 * A0 + 0.07667 * A0 * A0 - 0.0736 * A0 * A0 *A0 + 0.01008* A0 * A0 * A0 * A0;
	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}

double BlackOil::muSatDindoruk_Christman() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	double tempMuDeadOil = PaS_to_cp(muDeadOil);
	double a1 = 1.0; double a2 = 4.740729*1e-4; double a3 = -1.023451*1e-2; double a4 = 6.600358*1e-1;
	double a5 = 1.075080*1e-3;  double a6 = 1.0; double a7 = -2.191172*1e-5; double a8 = -1.660981*1e-2;
	double a9 = 4.233179*1e-1; double a10 = -2.273945*1e-4;
	double A1 = a1 / exp(a2*tempRs) + a3 * pow(tempRs, a4) / exp(a5 * tempRs);
	double A2 = a6 / exp(a7*tempRs) + a8 * pow(tempRs, a9) / exp(a10 * tempRs);
	muSatOil = A1 * pow(tempMuDeadOil, A2);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}

double BlackOil::muSatOilFunc() {
	if (muSatOilMethod == "Standing") {
		return muSatStanding();
	}
	else if (muSatOilMethod == "Beggs-Robinson") {
		return muSatBeggs_Robinson();
	}
	else if (muSatOilMethod == "Bergman") {
		return muSatBergman();
	}
	else if (muSatOilMethod == "Aziz") {
		return muSatAziz();
	}
	else if (muSatOilMethod == "AlKhafaji") {
		return muSatAlKhafaji();
	}
	else if (muSatOilMethod == "Dindoruk_Christman") {
		return muSatDindoruk_Christman();
	}
	else {
		cout << "*************" << endl;
	}
}


double BlackOil::muUnderSatFunc() {
	muSatOilFunc();
	double C1 = 2.6; double C2 = 1.187;
	double C3 = -11.513; double C4 = -8.98 * 1e-5;
	double B = C1 * pow(P, C2) * exp(C3 + C4 * P);
	muUnsatOil = muSatOil * pow(P / Pb, B);
	return muUnsatOil;
}


double BlackOil::BgCalculation() {
	double Bg = (Psc / Tsc) * (T / P) * Z;
	return Bg;
}

double BlackOil::BgCalculation(double P, double T, double Z) {
	double Bg = (Psc / Tsc) * (T / P) * Z;
	return Bg;
}


double BlackOil::gammaOilCal() {
	gammaOil = 141.5 / (API + 131.5);
	return gammaOil;
}

double BlackOil::gammaOilCal(double API) {
	gammaOil = 141.5 / (API + 131.5);
	return gammaOil;
}


double BlackOil::stockTankMw(string corr) {

	if (corr == "PVTsim") {
		if (API <= 40.0) {
			MwOil = 630.0 - 10 * API;	// gr/mole or 1b/1bmole
		}
		else {
			MwOil = 73110.0 * pow(API, -1.562);
		}
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate Stovk tank oil molecular weight.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
	return MwOil;
}

void BlackOil::Pc_Tc_Sutton() {
	gasTpc = R_to_K(169.2 + 349.5 * gammaGas - 74.0*gammaGas * gammaGas);
	gasPpc = psi_to_Pas(756.8 - 131.0 * gammaGas - 3.6 * gammaGas * gammaGas);
}

void BlackOil::gasCriticalProps(string corr) {
	if (corr == "Sutton") {
		Pc_Tc_Sutton();
	}
	else {
		string msg;
		msg = "Undefined function: " + corr + " to caculate gas specific gravity.";
		cout << "E0001	" << msg << endl;
		exit(1);
	}
}

double BlackOil::muGasDempsey() {
	double tempT = K_to_F(T);
	double Tpr = (tempT + 495.67) / (175.59 + 307.97 * gammaGas);
	double tempP = Pas_to_psi(P);
	double Ppr = tempP / (700.55 - 47.44 * gammaGas);
	double a0 = -2.46211820; double a1 = 2.97054714; double a2 = 2.86264054*1e-1; double a3 = 8.05420533*1e-3;
	double a4 = 2.80860949; double a5 = -3.49803305; double a6 = 3.60373020*1e-1; double a7 = -1.04432413*1e-2;
	double a8 = -7.93385684*1e-1; double a9 = 1.39643306; double a10 = -1.49144925*1e-1; double a11 = 4.41015512*1e-3;
	double a12 = 8.39387178*1e-2; double a13 = -1.86408848*1e-1; double a14 = 2.03367881*1e-2; double a15 = -6.09579263*1e-4;

	double temp1 = a0 + a1 * Ppr + a2 * Ppr * Ppr + a3 * Ppr * Ppr * Ppr +
		Tpr * (a4 + a5 * Ppr + a6 * Ppr * Ppr + a7 * Ppr * Ppr * Ppr);
	double temp2 = Tpr * Tpr * (a8 + a9 * Ppr + a10 * Ppr * Ppr + a11 * Ppr * Ppr * Ppr) +
		Tpr * Tpr * Tpr * (a12 + a13 * Ppr + a14 * Ppr * Ppr + a15 * Ppr * Ppr * Ppr);
	double lnf = temp1 + temp2;

	double b0 = 1.11231913*1e-2; double b1 = 1.67726604*1e-5; double b2 = 2.11360496 * 1e-9; double b3 = -1.09485050*1e-4;
	double b4 = -6.40316395*1e-8; double b5 = -8.99374533*1e-11; double b6 = 4.57735189*1e-7; double b7 = 2.12903390*1e-10;
	double b8 = 3.97732249*1e-13;

	double temp3 = b0 + b1 * tempT + b2 * tempT * tempT + b3 * MwGas + b4 * tempT * MwGas + b5 * tempT * tempT * MwGas
		+ b6 * MwGas * MwGas + b7 * tempT * MwGas * MwGas + b8 * tempT * tempT * MwGas * MwGas;

	muGas = (temp3 / Tpr) * exp(lnf);
	muGas = cp_to_PaS(muGas);
	return muGas;
}

double BlackOil::muGasLee_Gunzalez() {
	double temp_MwGas = MwGas * 1000.0;
	double temp_dGas = dGas / 1000.0;
	double tempT = K_to_R(T);
	double A1 = (9.379 + 0.01607 * temp_MwGas) * pow(tempT, 1.5) / (209.2 + 19.26 * temp_MwGas + tempT);
	double A2 = 3.448 + 986.6 / tempT + 0.01009 * temp_MwGas;
	double A3 = 2.447 - 0.2224 * A2;
	muGas = (A1 * 1e-4) * exp(A2 * pow(temp_dGas, A3));
	muGas = cp_to_PaS(muGas);
	return muGas;
}

double BlackOil::muGasLee_Gunzalez(double T, double MwGas, double dGas) {
	double tempT = K_to_R(T);
	double A1 = (9.379 + 0.01607 * MwGas) * pow(tempT, 1.5) / (209.2 + 19.26 * MwGas + tempT);
	double A2 = 3.448 + 986.6 / tempT + 0.01009 * MwGas;
	double A3 = 2.447 - 0.2224 * A2;
	muGas = (A1 * 1e-4) * exp(A2 * pow(dGas, A3));
	muGas = cp_to_PaS(muGas);
	return muGas;
}



double BlackOil::muGasLucas() {
	// 1 < Tpr < 40 & 0 < Ppr < 100
	double tempPpc = Pas_to_psi(gasPpc);
	double tempTpc = K_to_R(gasTpc);
	double tempP = Pas_to_psi(P);
	double tempT = K_to_R(T);
	double Tpr = tempT / tempTpc;
	double Ppr = tempP / tempPpc;
	double A1 = (1.245*1e-3) * exp(5.1726 * pow(Tpr, -0.3286)) / Tpr;
	double A2 = A1 * (1.6553 * Tpr - 1.2723);
	double A3 = 0.4489 * exp(3.0578 * pow(Tpr, -37.7332)) / Tpr;
	double A4 = 1.7368 * exp(2.2310 * pow(Tpr, -7.6351)) / Tpr;
	double A5 = 0.9425 * exp(-0.1853 * pow(Tpr, 0.4489));
	double sai = 9490.0 * pow(tempTpc / (MwGas * MwGas * MwGas * tempPpc * tempPpc * tempPpc * tempPpc), 1.0 / 6.0);
	double eta_gsc = (0.807 * pow(Tpr, 0.618) - 0.357 * exp(-0.449 * Tpr) + 0.340 * exp(-4.058 * Tpr) + 0.018) / sai;
	muGas = eta_gsc * (1 + A1 * pow(Ppr, 1.3088) / (A2 * pow(Ppr, A5) + 1.0 / (1 + A3 * pow(Ppr, A4))));
	muGas = cp_to_PaS(muGas);
	return muGas;
}

double BlackOil::muGasLucas(double gasPpc, double gasTpc, double P, double T, double MwGas) {
	// 1 < Tpr < 40 & 0 < Ppr < 100
	double tempPpc = Pas_to_psi(gasPpc);
	double tempTpc = K_to_R(gasTpc);
	double tempP = Pas_to_psi(P);
	double tempT = K_to_R(T);
	double Tpr = tempT / tempTpc;
	double Ppr = tempP / tempPpc;
	double A1 = (1.245*1e-3) * exp(5.1726 * pow(Tpr, -0.3286)) / Tpr;
	double A2 = A1 * (1.6553 * Tpr - 1.2723);
	double A3 = 0.4489 * exp(3.0578 * pow(Tpr, -37.7332)) / Tpr;
	double A4 = 1.7368 * exp(2.2310 * pow(Tpr, -7.6351)) / Tpr;
	double A5 = 0.9425 * exp(-0.1853 * pow(Tpr, 0.4489));
	double sai = 9490.0 * pow(tempTpc / (MwGas * MwGas * MwGas * tempPpc * tempPpc * tempPpc * tempPpc), 1.0 / 6.0);
	double eta_gsc = (0.807 * pow(Tpr, 0.618) - 0.357 * exp(-0.449 * Tpr) + 0.340 * exp(-4.058 * Tpr) + 0.018) / sai;
	muGas = eta_gsc * (1 + A1 * pow(Ppr, 1.3088) / (A2 * pow(Ppr, A5) + 1.0 / (1 + A3 * pow(Ppr, A4))));
	muGas = cp_to_PaS(muGas);
	return muGas;
}

double BlackOil::muGasFunc() {
	if (muGasMethod == "Lucas") {
		return muGasLucas();
	}
	else if (muGasMethod == "Lee-Gunzalez") {
		return muGasLee_Gunzalez();
	}
	else if (muGasMethod == "Dempsey") {
		return muGasDempsey();
	}
}



double BlackOil::ZFactor_PVTsim() {
	double tempPpc = Pas_to_psi(gasPpc);
	double tempTpc = K_to_R(gasTpc);
	double tempP = Pas_to_psi(P);
	double tempT = K_to_R(T);
	double Tpr = tempT / tempTpc;
	double Ppr = tempP / tempPpc;
	double t = 1 / Tpr;

	double fy = 1.0;
	double y = 0.001;
	double dfy;
	double ynew;
	while (abs(fy) > tol1) {

		fy = -0.06125 * Ppr * t * exp(-1.2 * (1 - t) * (1 - t)) +
			((y + y * y + y * y * y - y * y * y * y) / ((1 - y) * (1 - y) * (1 - y))) -
			(14.76 * t - 9.76 * t * t + 4.58 * t * t * t) * y * y +
			(90.7 * t - 242.2 * t * t + 42.4 * t * t * t) * pow(y, 2.18 + 2.82 * t);
		dfy = (1 + 4 * y + 4 * y * y - 4 * y * y * y + y * y * y * y) / ((1 - y) * (1 - y) * (1 - y) * (1 - y)) -
			(29.52 * t - 19.52 * t * t + 9.16 * t * t * t) * y + (2.18 + 2.82 * t)*(90.7 * t - 242.2 * t * t + 42.4 * t * t * t) *
			pow(y, 1.18 + 2.82 * t);

		ynew = y - fy / dfy;
		y = ynew;
	}

	Z = 0.06125 * Ppr * t * exp(-1.2 * (1 - t) * (1 - t)) / y;

	return Z;
}

double BlackOil::ZFactor_PVTsim(double P, double T, double gasPpc, double gasTpc) {
	double tempPpc = Pas_to_psi(gasPpc);
	double tempTpc = K_to_R(gasTpc);
	double tempP = Pas_to_psi(P);
	double tempT = K_to_R(T);
	double Tpr = tempT / tempTpc;
	double Ppr = tempP / tempPpc;
	double t = 1 / Tpr;

	double fy = 1.0;
	double y = 0.001;
	double dfy;
	double ynew;
	while (abs(fy) > tol1) {

		fy = -0.06125 * Ppr * t * exp(-1.2 * (1 - t) * (1 - t)) +
			((y + y * y + y * y * y - y * y * y * y) / ((1 - y) * (1 - y) * (1 - y))) -
			(14.76 * t - 9.76 * t * t + 4.58 * t * t * t) * y * y +
			(90.7 * t - 242.2 * t * t + 42.4 * t * t * t) * pow(y, 2.18 + 2.82 * t);
		dfy = (1 + 4 * y + 4 * y * y - 4 * y * y * y + y * y * y * y) / ((1 - y) * (1 - y) * (1 - y) * (1 - y)) -
			(29.52 * t - 19.52 * t * t + 9.16 * t * t * t) * y + (2.18 + 2.82 * t)*(90.7 * t - 242.2 * t * t + 42.4 * t * t * t) *
			pow(y, 1.18 + 2.82 * t);

		ynew = y - fy / dfy;
		y = ynew;
	}

	Z = 0.06125 * Ppr * t * exp(-1.2 * (1 - t) * (1 - t)) / y;

	return Z;
}


double BlackOil::gasMolecularWeight() {
	MwGas = gammaGas * 28.97;;
	return MwGas;
}

double BlackOil::GasDensity() {
	dGas = P * MwGas / (Z * Rgas * T);
	return dGas;
}

double BlackOil::GasDensity(double P, double T, double Z) {
	dGas = P * MwGas / (Z * Rgas * T);
	return dGas;
}

double BlackOil::oilDensityBubble(double gammaOil, double Rsb, double gammaGas, double Bob) {
	dOil = (62.4 * gammaOil + 0.01357 * Rsb * gammaGas) / Bob;
	return dOil;
}

double BlackOil::oilDensityBubble() {
	double tempRs = scmPerscm_To_scfPerSTB(Rs);
	dOil = ((62.4 * gammaOil + 0.01357 * tempRs * gammaGas) / Bob) / 62.4 * 1000.0;
	return dOil;
}

double BlackOil::oilDensityAboveBubble() {
	oilDensityBubble();
	dOilAboveBubble = dOil * exp(Co * (P - Pb));
	return dOilAboveBubble;
}

double BlackOil::muWaterFunc() {
	double tempT = K_to_F(T);
	double tempP = Pas_to_psi(P);
	double A0 = 109.574;	double A1 = -8.40564;		double A2 = 0.313314;			double A3 = 8.72213 * 1e-3;
	double B0 = -1.12166;	double B1 = 2.63951 * 1e-2; double B2 = -6.79461 * 1e-4;	double B3 = -5.47119 * 1e-5;
	double B4 = 1.55586 * 1e-6;
	double A = A0 + A1 * waterSal + A2 * waterSal * waterSal + A3 * waterSal * waterSal * waterSal;
	double B = B0 + B1 * waterSal + B2 * waterSal * waterSal + B3 * waterSal * waterSal * waterSal +
		B4 * waterSal * waterSal * waterSal * waterSal;
	double muW1 = A * pow(tempT, B);
	muSatOil = muW1 * (0.9994 + 4.0295 * 1e-5 * tempP + 3.1062 * 1e-9 * tempP * tempP);
	muSatOil = cp_to_PaS(muSatOil);
	return muSatOil;
}




void BlackOil::set_T(double var) {
	T = var;
}
void BlackOil::set_P(double var) {
	P = var;
}
void BlackOil::set_Rs(double var) {
	Rs = var;
}
void BlackOil::set_Bob(double var) {
	Bob = var;
}
void BlackOil::set_Bg(double var) {
	Bg = var;
}
void BlackOil::set_gammaGas(double var) {
	gammaGas = var;
}
void BlackOil::set_gammaOil(double var) {
	gammaOil = var;
}
void BlackOil::set_API(double var) {
	API = var;
}
void BlackOil::set_Pb(double var) {
	Pb = var;
}
void BlackOil::set_MwOil(double var) {
	MwOil = var;
}
void BlackOil::set_MwGas(double var) {
	MwGas = var;
}
void BlackOil::set_Z(double var) {
	Z = var;
}
void BlackOil::set_gasTpc(double var) {
	gasTpc = var;
}
void BlackOil::set_gasPpc(double var) {
	gasPpc = var;
}
void BlackOil::set_muDeadOil(double var) {
	muDeadOil = var;
}
void BlackOil::set_muSatOil(double var) {
	muSatOil = var;
}
void BlackOil::set_muGas(double var) {
	muGas = var;
}
void BlackOil::set_muOil(double var) {
	muOil = var;
}
void BlackOil::set_dGas(double var) {
	dGas = var;
}
void BlackOil::set_dOil(double var) {
	dOil = var;
}



double BlackOil::get_T() {
	return T;
}
double BlackOil::get_P() {
	return P;
}
double BlackOil::get_Rs() {
	return Rs;
}
double BlackOil::get_Bob() {
	return Bob;
}
double BlackOil::get_Bg() {
	return Bg;
}
double BlackOil::get_gammaGas() {
	return gammaGas;
}
double BlackOil::get_gammaOil() {
	return gammaOil;
}
double BlackOil::get_API() {
	return API;
}
double BlackOil::get_Pb() {
	return Pb;
}
double BlackOil::get_MwOil() {
	return MwOil;
}
double BlackOil::get_MwGas() {
	return MwGas;
}
double BlackOil::get_Z() {
	return Z;
}
double BlackOil::get_gasTpc() {
	return gasTpc;
}
double BlackOil::get_gasPpc() {
	return gasPpc;
}
double BlackOil::get_muDeadOil() {
	return muDeadOil;
}
double BlackOil::get_muSatOil() {
	return muSatOil;
}
double BlackOil::get_muGas() {
	return  muGas;
}
double BlackOil::get_muOil() {
	return muOil;
}
double BlackOil::get_dGas() {
	return dGas;
}
double BlackOil::get_dOil() {
	return dOil;
}

