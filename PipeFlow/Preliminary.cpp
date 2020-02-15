#include "Preliminary.h"



Preliminary::Preliminary()
{
}


Preliminary::~Preliminary()
{
}


double Preliminary::superficialVelocity(double q, double A) {
	return q / A;
}

double Preliminary::insituVelocity(double q, double A_phase) {
	return q / A_phase;
}

double Preliminary::twoPhaseFlow(double v_sL, double v_sg) {
	return v_sL + v_sg;
}

double Preliminary::volumeFractionFlowRate(double s, double m) {
	return s / m;
}

double Preliminary::massFractionFlowRate(double rho_1, double q_1, double rho_2, double q_2) {
	return (rho_1 * q_1 / (rho_1 * q_1 + rho_2 * q_2));
}

double Preliminary::massFlux(double w_1, double w_2, double A) {
	return (w_1 + w_2) / A;
}

double Preliminary::phaseMassFlux(double x, double Gm) {
	return x * Gm;
}

double Preliminary::slipVelocity(double v_g, double v_L) {
	return v_g - v_L;
}

double Preliminary::insituVelocityIBynsituVolumeFraction(double V_s, double f) {
	return V_s / f;
}

/*double Preliminary::frictionalPressureGradient(double f, double v, double rho, double D) {
	double dPdH_F = f * v * v * rho / (2 * D);
	return dPdH_F;
}*/

/*double Preliminary::HydrostaticPressureGradient(double rho, double theta) {
	double dPdH_H = rho * g_acc * sin(theta);
	return dPdH_H;
}*/

/*double Preliminary::accelerationPressureGradient(double rho, double v, double dvdH) {
	double dPdH_A = rho * v * dvdH;
	return dPdH_A;
}*/

/*double Preliminary::total_dPdH(double dPdH_H, double dPdH_F, double dPdH_A) {
	return dPdH_A + dPdH_F + dPdH_H;
}*/

double Preliminary::rhoMixtureVol(double rho_L, double rho_g, double f_L) {
	return f_L * rho_L + (1 - f_L) * rho_g;
}

double Preliminary::rhoMixtureMass(double rho_L, double rho_g, double x_g) {
	return rho_g * rho_L / (x_g * rho_L + (1 - x_g) * rho_g);
}

double Preliminary::muMixtureCcchitti(double mu_L, double mu_g, double x_g) {
	return x_g * mu_g + (1 - x_g) * mu_L;
}

double Preliminary::muMixtureDuckler(double mu_L, double mu_g, double x_g, double rho_L, double rho_g, double rho_m) {
	return rho_m * (x_g * mu_g / rho_g + (1 - x_g) * mu_L / rho_L);
}

double Preliminary::muMixtureMcAdams(double mu_L, double mu_g, double x_g) {
	return mu_g * mu_L / (x_g * mu_L + (1 - x_g) * mu_g);
}

double Preliminary::muMixtureVol(double mu_L, double mu_g, double f_L) {
	return f_L * mu_L + (1 - f_L) * mu_g;
}

double Preliminary::v_sgBubblyFlowAnsari(double v_sL, double v_inf) {
	double v_sg = 0.25 * v_inf + 0.333 * v_sL;
	return v_sg;
}

double Preliminary::Harmathy(double rho_L, double rho_g, double sigma) {
	double v_inf = 1.53 * sqrt(sqrt(g_acc * (rho_L - rho_g) * sigma / (rho_L * rho_L)));
	return v_inf;
}

double Preliminary::v_sgBubbleToSlugZuber(double v_sL, double v_inf) {
	double v_sg = 0.429 * v_sL + 0.357 * v_inf;
	return v_sg;
}

double Preliminary::v_sgBubbleToSlugTaitel(double v_sL, double v_inf) {
	double v_sg = 0.33 * v_sL + 0.25 * v_inf;
	return v_sg;
}

double Preliminary::riseVelocityTaylorBuble(double rho_L, double rho_g, double d) {
	double v_infT = 0.35 * sqrt(g_acc * d * (rho_L - rho_g) / rho_L);
	return v_infT;
}

bool Preliminary::bubbleFlowExistenceTaitel(double d, double rho_L, double rho_g, double sigma) {
	double RHS = 19.1 * sqrt(sigma / (g_acc * (rho_L - rho_g)));
	if (d < RHS) {
		return false;
	}
	else {
		return true;
	}
}

bool Preliminary::dispersedBubblyFlowTaitelSoham(double v_m, double v_sg, double f, double d, double rho_L, double rho_g, double sigma) {
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

bool Preliminary::transitionToChurn(double v_sg, double v_sL) {
	if (v_sg > 1.08 * v_sL) {
		return true;
	}
	else {
		return false;
	}
}

bool Preliminary::transitionToChurn(double v_sg, double rho_L, double rho_g, double sigma) {
	double v_sg_cal = 3.1 * sqrt(sqrt(g_acc * sigma * (rho_L - rho_g) / (rho_g * rho_g)));
	if (v_sg > v_sg_cal) {
		return v_sg_cal;
	}
	else {
		return false;
	}
}

double Preliminary::BubblyFlowInsituGasVelocity(double v_m, double v_d, double f_g) {
	double Co = 1.02;		// small pipe diameter for the large one ( > 127 mm) Co = 2.0
							// we only have deal with the small ones
	double v_g = Co * v_m + v_d / f_g;
	return v_g;
}

double Preliminary::BubblyFlowTerminalDriftVelocity(double v_inf, double f_g) {
	//double n;
	//n = 0.0;	// for vertical flow
	//double v_d = v_inf * f_g * pow(1 - f_g, n);
	double v_d = v_inf * f_g;
	return v_d;
}

double Preliminary::BubblyFlowF_g(double v_m, double v_inf, double v_sg, string model) {
	double Co = 1.2;		// small pipe diameter for the large one ( > 127 mm) Co = 2.0
	double df_g, f_g, f_g_old;		// we only have deal with the small ones
	f_g_old = 0.1;
	if ((model == "Ansari") || (model == "ansari")) {

		double err = 1.0;
		while (err > 1e-8) {

			f_g = v_sg / (Co * v_m + v_inf * sqrt(1 - f_g_old));
			df_g = v_sg / sqrt(1 - f_g_old) / (Co * v_m + v_inf * sqrt(1 - f_g_old)) / (Co * v_m + v_inf * sqrt(1 - f_g_old));
			f_g = f_g_old - f_g / df_g;
			err = abs((f_g - f_g_old) / f_g);
			f_g = f_g_old;
		}

	}
	else {
		double f_g = v_sg / (Co * v_m + v_inf);
	}
	return f_g;
}

double Preliminary::slugFlowV_sg(double v_sL, double v_inf) {
	double v_sg = 0.429 * v_sL + 0.357 * v_inf;
	return v_sg;
}

double Preliminary::slugFlowFallingLiqFilmVelocity(double d, double f_gTB) {
	double v_LTB = 9.916 * sqrt(g_acc * d * (1 - sqrt(f_gTB)));
	return v_LTB;
}

double Preliminary::slugFlow_f_gLS(double v_sg, double v_m) {
	double f_gLS = v_sg / (2.65 * v_m + 0.425);
	return f_gLS;
}

void Preliminary::BubbleFlowMixtureProps(double rho_L, double rho_g, double f_L, double& rho_m, double v_m, double d, double mu_m) {
	rho_m = rho_L * f_L + rho_g * (1 - f_L);
	double Re = ReynoldNumber(rho_m, v_m, d, mu_m);

}


double Preliminary::colebrookWhiteMoody(double e, double d, double Re) {
	double f_d;
	double f_d_old = 1.0 / ((2.0 * log10(e / (3.7 * d))) * (2.0 * log10(e / (3.7 * d))));
	double err = 1.0;
	while (err > 1e-8)
	{
		f_d = 1.0 / ((2.0 * log10(e / (3.7 * d) + 2.51 / Re / sqrt(f_d_old))) *
			(2.0 * log10(e / (3.7 * d) + 2.51 / Re / sqrt(f_d_old))));
		err = abs((f_d - f_d_old) / f_d_old);
		f_d_old = f_d;
	}
	return f_d;
}
double Preliminary::ReynoldNumber(double rho, double v, double d, double mu) {
	double Re = rho * v * d / mu;
	return Re;
}


double Preliminary::singlePhase_dPdH_PE(double rho, double theta) {
	double dpdz_pe = rho * g_acc * sin(theta);
	return dpdz_pe;
}
double Preliminary::singlePhase_dPdH_KE(double rho, double q, double D1, double D2) {
	double dpdz_ke = 8.0 * rho * q * q / piNum / piNum * ((1.0 / pow(D2, 4.0)) - (1.0 / pow(D1, 4.0)));
	return dpdz_ke;
}
double Preliminary::singlePhase_dPdH_F(double rho, double u, double D, double f) {
	double dpdz_f = f / 2.0 * rho * u * u / D;
	return dpdz_f;
}