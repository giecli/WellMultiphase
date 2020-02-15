#include "twoPhaseModels.h"



twoPhaseModels::twoPhaseModels()
{
}


twoPhaseModels::~twoPhaseModels()
{
}

bool twoPhaseModels::bubbleFlow(double rho_L, double rho_g, double sigma, double d, double f_g) {
	bool bubbleFlowRegime;
	double d_min = 19.1 * sqrt(sigma / (g_acc * (rho_L - rho_g)));
	/*if (d < d_min) {
		bubbleFlowRegime = false;
		return bubbleFlowRegime;
	}
	else*/ if (f_g < 0.25) {
		bubbleFlowRegime = true;
		return bubbleFlowRegime;
	}
	else if (f_g >= 0.25) {
		bubbleFlowRegime = false;
		return bubbleFlowRegime;
	}
}

bool twoPhaseModels::dispersedBubbleFlow(double v_m, double v_sg, double f, double d, double rho_L, double rho_g, double sigma, double f_g) {
	// f liquid friction factor
	if (f_g > 0.52) {
		return false;
	}
	else {
		double err = 1.0;
		double ff, dff;
		double v_m_cal_new, v_m_cal_old;
		v_m_cal_new = v_m;
		while (err > 1e-8) {
			ff = 2.0 * pow(v_m_cal_new, 1.2) * pow(f / (2.0 * d), 0.4) * pow(rho_L / sigma, 0.6) *
				sqrt(0.4 * sigma / (g_acc * (rho_L - rho_g))) - 0.725 - 4.15 * sqrt(v_sg / v_m_cal_new);
			dff = 2.4 * pow(v_m_cal_new, 1.2) * pow(f / (2.0 * d), 0.4) * pow(rho_L / sigma, 0.6) *
				sqrt(0.4 * sigma / (g_acc * (rho_L - rho_g))) - 4.15 *
				sqrt(v_sg) * (-1.0 / 2.0) / sqrt(v_m_cal_new * v_m_cal_new * v_m_cal_new);
			v_m_cal_old = v_m_cal_new;
			err = abs((v_m_cal_old - v_m_cal_new) / v_m_cal_old);
		}
		if (v_m > v_m_cal_new) {
			return true;
		}
		else {
			return false;
		}
	}
}

bool twoPhaseModels::transitionToAnnularFlow(double rho_L, double rho_g, double sigma, double v_sg) {
	double v_inf = 1.53 * sqrt(sqrt(g_acc * (rho_L - rho_g) * sigma / (rho_L * rho_L)));
	if (v_sg > v_inf) {
		return true;
	}
	else {
		return false;
	}
}

bool twoPhaseModels::transitionToChurnFlow(double f_g) {
	if (f_g > 0.76) {
		return true;
	}
	else {
		return false;
	}
}

string twoPhaseModels::ansariPattern(double rho_L, double rho_g, double sigma, double d, double f_g, double v_m, double v_sg, double f) {
	string pattern = "Non";

	if (bubbleFlow(rho_L, rho_g, sigma, d, f_g)) {
		pattern = "bubble flow";
	}
	else {
		if (dispersedBubbleFlow(v_m, v_sg, f, d, rho_L, rho_g, sigma, f_g)) {
			pattern = "dispersed bubble flow";
		}
		else {
			if (transitionToChurnFlow(f_g)) {
				if (transitionToAnnularFlow(rho_L, rho_g, sigma, v_sg)) {
					pattern = "annular flow";
				}
				else {
					pattern = "churn flow";
				}
			}
			else {
				pattern = "slug flow";
			}
		}
	}
	return pattern;
}


double twoPhaseModels::dpdzEleBubbleFlow(double rho_m, double theta) {
	double dpdz = rho_m * g_acc * sin(theta);
	double aaa = sin(theta);
	return dpdz;
}
double twoPhaseModels::dpdzFricBubbleFlow(double f_TP, double rho_m, double v_m, double d) {
	double dzdp = f_TP * rho_m * v_m * v_m / (2.0 * d);
	return dzdp;
}




double twoPhaseModels::ansariPressureGradient(string pattern, double &f_L, double &f_TP, double &rho_m, double &mu_m,
	double &v_m, double d, double theta,
	double rho_L, double rho_g, double mu_L, double mu_g, double v_sL, double v_sg, double &dv_mdz, double e,
	double sigma) {
	double dpdz;
	if ((pattern == "bubble flow") || (pattern == "dispersed bubble flow")) {
		v_m = v_sg + v_sL;
		f_L = v_sL / v_m;
		rho_m = rhoMixtureVol(rho_L, rho_g, f_L);
		mu_m = muMixtureVol(mu_L, mu_g, f_L);
		double Re = ReynoldNumber(rho_m, v_m, d, mu_m);
		f_TP = colebrookWhiteMoody(e, d, Re);
		double dpdzEle = dpdzEleBubbleFlow(rho_m, theta);
		double dpdzFri = dpdzFricBubbleFlow(f_TP, rho_m, v_m, d);
		dpdz = -dpdzEle - dpdzFri;
	}
	else if (pattern == "slug flow") {
		v_m = v_sg + v_sL;
		double v_inf = Harmathy(rho_L, rho_g, sigma);
		double v_infT = 0.35 * sqrt(g_acc * d * (rho_L - rho_g) / rho_L);
		double f_gLS = v_sg / (2.65 * v_m + 0.425);
		double f_LLS = 1.0 - f_gLS;
		double v_gLS = 1.2 * v_m + v_inf * sqrt(f_LLS);
		double v_TB = 1.2 * v_m + v_infT;
		double A = v_TB * f_gLS + (1.0 - f_gLS) * (v_m - f_gLS * v_inf * sqrt(1.0 - f_gLS));
		double err = 1.0;
		double f_gTB_old = 2.0 * f_gLS;
		double v_LTB;
		double f_LTB, f_gTB;
		while (err > 1e-8) {
			v_LTB = 9.916 * sqrt(g_acc * d * (1.0 - sqrt(f_gTB_old)));
			f_LTB = (v_TB - A) / (v_LTB + v_TB);
			f_gTB = 1.0 - f_LTB;
			err = abs((f_gTB - f_gTB_old) / f_gTB_old);
			f_gTB_old = f_gTB;
		}
		double v_LLS = v_TB - (v_TB + v_LTB) * f_LTB / f_LLS;
		double v_gTB = v_TB - (v_TB - v_gLS) * f_gLS / (1.0 - f_LLS);
		double B = (v_sg - v_gLS * (1.0 - f_LLS)) / (v_gTB * (1.0 - f_LTB) - v_gLS * (1.0 - f_LLS));
		double rho_LS = rho_L * f_LLS + rho_g * (1.0 - f_LLS);
		double mu_LS = mu_L * f_LLS + mu_g * (1.0 - f_LLS);
		//rho_m = rhoMixtureVol(rho_L, rho_g, f_L);
		//mu_m = muMixtureVol(mu_L, mu_g, f_L);
		double dpdzEle = dpdzEleBubbleFlow((1.0 - B) * rho_LS + B * rho_g, theta);
		double Re = ReynoldNumber(rho_LS, v_m, d, mu_LS);
		double fri_LS = colebrookWhiteMoody(e, d, Re);
		double dpdzFri = dpdzFricBubbleFlow(fri_LS, rho_LS * (1.0 - B), v_m, d);
		dpdz = -dpdzEle - dpdzFri;
	}

	else if (pattern == "annular flow") {

		double v_sgc = 10000.0 * v_sg * mu_g / sigma * sqrt(rho_g / rho_L);
		double E = 1.0 - exp(-0.125 * (v_sgc - 1.5));
		double lambda_LC = 1.0 - v_sg / (v_sg + E * v_sL);
		double rho_c = rho_L * lambda_LC + rho_g * (1.0 - lambda_LC);
		double mu_c = mu_L * lambda_LC + mu_g * (1.0 - lambda_LC);
		double v_sc = v_sg + E * v_sL;

		double Z_cor, v_c, f_c, Re_c;
		double err = 1.0;
		double delta_old = 0.0001;
		double delta;
		double f;

		double Re_LF, f_LF;

		double fDelta, dfDelta;
		double a1, a2;
		while (err > 1e-8) {
			v_c = v_sc / ((1.0 - 2.0 * delta_old) * (1.0 - 2.0 * delta_old));
			Re_c = ReynoldNumber(rho_c, v_c, 1 - 2.0 * delta_old, mu_c);
			//f = colebrookWhiteMoody(e, 1.0 - 2.0 * delta_old, Re_c);
			f = 0.32 / sqrt(sqrt(Re_c));
			if (E >= 0.9) {
				f_c = f * (1.0 + 300.0 * delta_old);
			}
			else {
				f_c = f * (1.0 + 24.0 * cbrt(rho_L / rho_g) * delta_old);
			}
			Re_LF = ReynoldNumber(rho_L, v_sL, d, mu_L);
			f_LF = 4.0 * colebrookWhiteMoody(e, d, Re_LF);
			a1 = f_c * rho_c * v_c * v_c / (2.0 * d * delta_old * (1.0 - delta_old) * (1.0 - 2.0 * delta_old));
			a2 = 2.0 * f_LF * v_sL * v_sL * (1.0 - E) * (1.0 - E) * rho_L / (64.0 * d * delta_old * delta_old * delta_old *
				(1.0 - delta_old) * (1.0 - delta_old) * (1.0 - delta_old));
			fDelta = a1 - g_acc * (rho_L - rho_c) - a2;
			dfDelta = -a1 * ((1.0 - delta_old) * (1.0 - 2.0 * delta_old) - delta_old * (1.0 - 2.0 * delta_old) - 2.0 * delta_old * (1.0 - delta_old)) /
				(delta_old * (1.0 - delta_old) * (1.0 - 2.0 * delta_old)) + a2 * (3.0 * delta_old * delta_old *
				(1.0 - 2.0 * delta_old) * (1.0 - 2.0 * delta_old) * (1.0 - 2.0 * delta_old) - 6 * delta_old * delta_old * delta_old *
					(1.0 - 2.0 * delta_old) * (1.0 - 2.0 * delta_old)) /
					(delta_old * delta_old * delta_old *
				(1.0 - delta_old) * (1.0 - delta_old) * (1.0 - delta_old));
			delta = delta_old - fDelta / dfDelta;
		}
		double dpdzFri, dpdzEle;
		dpdzFri = f_c * v_c * v_c * rho_c / (2.0 * d * (1.0 - 2.0 * delta));
		dpdzEle = g_acc * rho_c;
		dpdz = -dpdzEle - dpdzFri;

	}
	return dpdz;
}

void twoPhaseModels::FlowProfile(double Qg, double Ql, double Pw, int gridNum, double T, double P_boundary, double gammaGas, double API,
	double Rs, double Co, double e, double d, double dL, double theta, string units, string Pseudocritical_method, string MwMethod,
	string muGasMethod, string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod, string flowDir, 
	double* PAns, double* L) {
	register int i;
	double A = piNum / 4.0 * d * d;
	PAns[0] = P_boundary;
	double P;
	double sign_dP;
	if (flowDir == "downstream") {
		sign_dP = 1.0;
		L[0] = dL * gridNum;
	}
	else {
		sign_dP = -1.0;
		L[0] = 0.0;
	}
	BlackOil blackOilObj(gammaGas, API, Rs, Co, units, Pseudocritical_method, MwMethod, muGasMethod,
		muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod);

	double ele = 0.0;
	for (i = 0; i < gridNum; i++) {
		L[i + 1] = L[i] - sign_dP * dL;
		double err = 1.0;
		PAns[i + 1] = PAns[i];
		while (err > 1e-8) {
			P = (PAns[i] + PAns[i + 1]) / 2.0;
			blackOilObj.set_T(T);
			blackOilObj.set_P(P);
			double rho_L, rho_g, sigma, f_g, v_m, v_sg, f, mu_L, mu_g;
			double Rsb = blackOilObj.RsStanding();
			blackOilObj.set_Rs(Rsb);
			double Qgb = (Rs - Rsb) * Ql;
			double dpdz;
			if (Qgb > 0) {
				double Bo = blackOilObj.BobFunc();
				rho_L = blackOilObj.oilDensityBubble();
				blackOilObj.muDeadOilFunc();
				mu_L = blackOilObj.SatOilViscosity();
				blackOilObj.ZFactor();
				rho_g = blackOilObj.GasDensity();
				mu_g = blackOilObj.muGasFunc();
				double Bg = blackOilObj.BgCalculation();
				v_sg = Bg * Qgb / A / Day_to_s(1.0);
				double v_sl = Ql * Bo / A / Day_to_s(1.0);
				v_m = v_sg + v_sl;
				sigma = dynePerCen_to_NPerM(30.0);
				double Re = ReynoldNumber(rho_L, v_sl, d, mu_L);
				f = colebrookWhiteMoody(e, d, Re);
				f_g = v_sg / v_m;
				string pattern = ansariPattern(rho_L, rho_g, sigma, d, f_g, v_m, v_sg, f);
				double f_L, f_TP, rho_m, mu_m, dv_mdz;
				dpdz = ansariPressureGradient(pattern, f_L, f_TP, rho_m, mu_m, v_m, d, theta,
					rho_L, rho_g, mu_L, mu_g, v_sl, v_sg, dv_mdz, e, sigma);
			}
			else {
				double Bo = blackOilObj.BoUnderSat();
				double dOil = blackOilObj.oilDensityAboveBubble();
				double v_sl = Ql * Bo / A / Day_to_s(1.0);
				double dpdz_pe = singlePhase_dPdH_PE(dOil, theta);
				double muOilAboveBunnle = blackOilObj.muUnderSatFunc();
				double dpdz_f = singlePhase_dPdH_F(dOil, v_sl, d, muOilAboveBunnle);
				double dpdz_ke = singlePhase_dPdH_KE(dOil, Ql, d, d);
			}
			err = abs(PAns[i] + sign_dP * dpdz * dL - PAns[i + 1]) / PAns[i + 1];
			PAns[i + 1] = sign_dP * dpdz * dL + PAns[i];
		}
	}
}



void twoPhaseModels::FlowProfile(double Qg, double Qo, double Qw, double waterSal, double Pw, int gridNum, double T, double P_boundary, double gammaGas, double API,
	double Rs, double Co, double e, double d, double dL, double theta, string units, string Pseudocritical_method, string MwMethod,
	string muGasMethod, string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod, string flowDir,
	double* PAns, double* L) {
	register int i;
	double A = piNum / 4.0 * d * d;
	/*double* PAns = new double[gridNum + 1];
	double* L = new double[gridNum + 1];*/
	for (i = 0; i < gridNum + 1; i++) {
		PAns[i] = 0.0;
	}
	PAns[0] = P_boundary;
	double P;
	double sign_dP;
	if (flowDir == "downstream") {
		sign_dP = 1.0;
		L[0] = dL * gridNum;
	}
	else {
		sign_dP = -1.0;
		L[0] = 0.0;
	}
	BlackOil blackOilObj(gammaGas, API, Rs, Co, waterSal, units, Pseudocritical_method, MwMethod, muGasMethod,
		muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod);

	double ele = 0.0;
	for (i = 0; i < gridNum; i++) {
		ele += dL;
		L[i + 1] = L[i] - sign_dP * dL;
		double err = 1.0;
		PAns[i + 1] = PAns[i];
		while (err > 1e-8) {
			P = (PAns[i] + PAns[i + 1]) / 2.0;
			blackOilObj.set_T(T);
			blackOilObj.set_P(P);
			double rho_L, rho_o, rho_W, rho_g, sigma, f_g, v_m, v_sg, f, mu_L, mu_o, mu_W, mu_g;
			double Rsb = blackOilObj.RsStanding();
			blackOilObj.set_Rs(Rsb);
			double Qgb = (Rs - Rsb) * Qo;
			double dpdz;
			if (Qgb > 0) {
				double Bo = blackOilObj.BobFunc();
				double Bw = blackOilObj.BwFunc();
				rho_o = blackOilObj.oilDensityBubble();
				/*rho_W = blackOilObj.waterDensityFunc();*/ //water density must be checked
				rho_W = 1050;
				blackOilObj.set_Pb(P);
				rho_L = rhoMixtureMass(rho_o, rho_W, Qw / (Qw + Qo));
				blackOilObj.muDeadOilFunc();
				mu_o = blackOilObj.SatOilViscosity();
				mu_W = blackOilObj.muWaterFunc();
				mu_L = muMixtureCcchitti(mu_o, mu_W, Qw / (Qw + Qo));
				blackOilObj.ZFactor();
				rho_g = blackOilObj.GasDensity();
				mu_g = blackOilObj.muGasFunc();
				double Bg = blackOilObj.BgCalculation();
				v_sg = Bg * Qgb / A / Day_to_s(1.0);
				double v_so = Qo * Bo / A / Day_to_s(1.0);
				double v_sw = Qw * Bw / A / Day_to_s(1.0);
				double v_sl = v_so + v_sw;
				v_m = v_sg + v_sl;
				sigma = dynePerCen_to_NPerM(31.6);
				double Re = ReynoldNumber(rho_L, v_sl, d, mu_L);
				f = colebrookWhiteMoody(e, d, Re);
				f_g = v_sg / v_m;
				string pattern = ansariPattern(rho_L, rho_g, sigma, d, f_g, v_m, v_sg, f);
				double f_L, f_TP, rho_m, mu_m, dv_mdz;
				dpdz = ansariPressureGradient(pattern, f_L, f_TP, rho_m, mu_m, v_m, d, theta,
					rho_L, rho_g, mu_L, mu_g, v_sl, v_sg, dv_mdz, e, sigma);
			}
			else {
				double Bo = blackOilObj.BoUnderSat();
				double Bw = blackOilObj.BwFunc();
				rho_o = blackOilObj.oilDensityAboveBubble();
				rho_W = blackOilObj.waterDensityFunc();
				rho_L = rhoMixtureMass(rho_o, rho_W, Qw / (Qw + Qo));
				double v_so = Qo * Bo / A / Day_to_s(1.0);
				double v_sw = Qw * Bw / A / Day_to_s(1.0);
				double v_sl = v_so + v_sw;
				double dpdz_pe = singlePhase_dPdH_PE(rho_L, theta);
				mu_o = blackOilObj.muUnderSatFunc();
				mu_W = blackOilObj.muWaterFunc();
				mu_L = muMixtureCcchitti(mu_o, mu_W, Qw / (Qw + Qo));
				double dpdz_f = singlePhase_dPdH_F(rho_L, v_sl, d, mu_L);
				double Ql = Qw + Qo;
				double dpdz_ke = singlePhase_dPdH_KE(rho_o, Ql, d, d);
				dpdz = -dpdz_f - dpdz_pe + dpdz_ke;
			}
			err = abs(PAns[i] + sign_dP * dpdz * dL - PAns[i + 1]) / PAns[i + 1];
			PAns[i + 1] = sign_dP * dpdz * dL + PAns[i];
		}
	}
	

}




