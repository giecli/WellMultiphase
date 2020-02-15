#pragma once
#include<string>

#include"global.h"

using namespace std;

class Preliminary
{
public:
	Preliminary();
	~Preliminary();
	// mixture properties
	double superficialVelocity(double q, double A);
	double insituVelocity(double q, double A_phase);
	double twoPhaseFlow(double v_sL, double v_sg);
	double volumeFractionFlowRate(double s, double m);
	double massFractionFlowRate(double rho_1, double q_1, double rho_2, double q_2);
	double massFlux(double w_1, double w_2, double A);
	double phaseMassFlux(double x, double Gm);
	double slipVelocity(double v_g, double v_L);
	double insituVelocityIBynsituVolumeFraction(double V_s, double f);
	double frictionalPressureGradient(double f, double v, double rho, double D);
	double HydrostaticPressureGradient(double rho, double theta);
	double HydrostaticPressureGradient(double rho, double v, double dvdH);
	double rhoMixtureVol(double rho_L, double rho_g, double f_L);
	double rhoMixtureMass(double rho_L, double rho_g, double x_g);
	double muMixtureCcchitti(double mu_L, double mu_g, double x_g);
	double muMixtureDuckler(double mu_L, double mu_g, double x_g, double rho_L, double rho_g, double rho_m);
	double muMixtureMcAdams(double mu_L, double mu_g, double x_g);
	double muMixtureVol(double mu_L, double mu_g, double f_L);

	// flow characteristic (f_g, f_L, V_sg, v_m, ....)
	double v_sgBubblyFlowAnsari(double v_sL, double v_inf);
	double Harmathy(double rho_L, double rho_g, double sigma);
	double v_sgBubbleToSlugZuber(double v_sL, double v_inf);
	double v_sgBubbleToSlugTaitel(double v_sL, double v_inf);
	double riseVelocityTaylorBuble(double rho_L, double rho_g, double d);
	bool bubbleFlowExistenceTaitel(double d, double rho_L, double rho_g, double sigma);
	bool dispersedBubblyFlowTaitelSoham(double v_m, double v_sg, double f, double d, double rho_L, double rho_g, double sigma);
	bool transitionToChurn(double v_sg, double v_sL);
	bool transitionToChurn(double v_sg, double rho_L, double rho_g, double sigma);
	double BubblyFlowInsituGasVelocity(double v_m, double v_d, double f_g);
	double BubblyFlowF_g(double v_m, double v_d, double f_g, string model);
	double BubblyFlowTerminalDriftVelocity(double v_inf, double f_g);

	double slugFlowV_sg(double v_sL, double v_inf);
	double slugFlowFallingLiqFilmVelocity(double d, double f_gTB);
	double slugFlow_f_gLS(double v_sg, double v_m);


	double singlePhase_dPdH_PE(double rho, double theta);
	double singlePhase_dPdH_KE(double rho, double q, double D1, double D2);
	double singlePhase_dPdH_F(double rho, double u, double D, double mu);
	//	double frictionalPressureGradient(double f, double v, double rho, double D);
	//	double HydrostaticPressureGradient(double rho, double theta);
	//	double accelerationPressureGradient(double rho, double v, double dvdH);

		// moody friction factor colebrook white
	double colebrookWhiteMoody(double e, double d, double Re);

	// calculation of flow regime propeties in pack
	double ReynoldNumber(double rho, double v, double d, double mu);
	void BubbleFlowMixtureProps(double rho_L, double rho_g, double f_L, double& rho_m, double v_m, double d, double mu);
private:
	string twoPhaseModel;
};

