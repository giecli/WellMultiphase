#pragma once

#include"Preliminary.h"
#include "BlackOil.h"
#include "Compositional.h"
#include "global.h"

class twoPhaseModels : Preliminary
{
public:
	twoPhaseModels();
	~twoPhaseModels();
	string ansariPattern(double rho_L, double rho_g, double sigma, double d, double f_g, double v_m, double v_sg, double f);
	double ansariPressureGradient(string pattern, double &f_L, double &f_TP, double &rho_m, double &mu_m,
		double &v_m, double d, double theta, double rho_L, double rho_g,
		double mu_L, double mu_g, double v_sL, double v_sg, double &dv_mdz, double e, double sigma);
	bool bubbleFlow(double rho_L, double rho_g, double sigma, double d, double f_g);
	bool dispersedBubbleFlow(double v_m, double v_sg, double f, double d, double rho_L, double rho_g, double sigma, double f_g);
	bool transitionToAnnularFlow(double rho_L, double rho_g, double sigma, double v_sg);
	bool transitionToChurnFlow(double f_g);

	double dpdzEleBubbleFlow(double rho_m, double theta);
	double dpdzFricBubbleFlow(double f_TP, double rho_m, double v_m, double d);
	void FlowProfile(double Qg, double Ql, double Pw, int gridNum, double T, double P_boundary, double gammaGas, double API,
		double Rs, double Co, double e, double d, double dL, double theta, string units, string Pseudocritical_method, string MwMethod,
		string muGasMethod, string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod, string dir,
		double* PAns, double* L);

	void FlowProfile(double Qg, double Qo, double Qw, double waterSal, double Pw, int gridNum, double T, double P_boundary, double gammaGas, double API,
		double Rs, double Co, double e, double d, double dL, double theta, string units, string Pseudocritical_method, string MwMethod,
		string muGasMethod, string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod, string dir,
		double* PAns, double* L);






};

