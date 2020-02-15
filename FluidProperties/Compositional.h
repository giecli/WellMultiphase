#pragma once
#include<string>
#include"global.h"
#include"mathFunc.h"
#include"BlackOil.h"
#include<Dense>


using namespace Eigen;
using namespace std;

class Compositional : public BlackOil
{
public:
	Compositional(double* Tc, double* Pc, double* ACF, double* ai, double* bi,
		int inNc, double T, double P, string inEOS);
	Compositional();
	~Compositional();
	double cubicEOS(double T, double P, double* comp, double* ai, double* bi,
		double& a, double& b, double* S, double& A, double& B, double** BIC, double& Z, string& phase);
	void RKEOS(double T, double P, double* Tc, double* Pc, double* a, double* b);
	void SRKEOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b);
	void GDEOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b);
	void PREOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b);
	void PR78EOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b);
	double fugacityCoef();
	void cubicFugacityCoef(double A, double B, double a, double b, double* bi,
		double* ai, double* Si, double Z, double* phi);
	void wilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* K);
	void dTwilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* dK);
	void dPwilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* dK);
	void mWilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* K);
	double cubicRootSelection(double ZA, double ZB, double B, double A, string& phase);

	void setNc(int N);
	void setEOS(string inEOS);

	int twoPhaseFlash_SS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv, double* x, double* y);
	int twoPhaseFlash_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv, double* x, double* y);
	void zeroInitializatinoCubicEOS(double a, double b, double* S, double A, double B, double Z, string phase);
	int twoPhaseFlash_QNN_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv, double* x, double* y);
	double twoPhaseStabilityFlash(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv);
	double twoPhaseStabilityFlash_BFGS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv);
	double twoPhaseStabilityFlash_BFGS_FullyImplicit(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv);
	int twoPhaseFlash_BFGS_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv, double* x, double* y);
	int twoPhaseFlash_BFGS_FullyImplicit(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double* K, double& nv, double* x, double* y);
	double enthalpyReal(double* x, double T, double P, double* a0, double* a1, double* a2, double* a3, double* a4,
		double Z, double A, double B, double phi);
	double enthalpyIdeal(double* x, double T, double* a0, double* a1, double* a2, double* a3, double* a4);
	double dZdTFunc(double Z, double A, double B, double T);
	double dphidTFunc(double dZdT, double A, double B, double Z, double T);
	double HeatCapacityIdeal(double* x, double T, double* a0, double* a1, double* a2, double* a3, double* a4);
	double enthalpyResedual(double T, double phi, double dphidT);
	double viscosity_LBC(double P, double T, double rho, double* x, double* Tc, double* Pc, double* Vc, double* Mw);
	double density_f(double P, double T, double MwMix, double Z);
	double molecularWight(double* x, double* Mw);
	double GasOilIFT(double rhoL, double rhoG, double* Par, double* x, double* y);
	int twoPhaseFlashEnv_BFGS_MSS_bubble(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF);
	int twoPhaseFlashEnv_BFGS_MSS_bubble2(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF);
	int twoPhaseFlashEnv_BFGS_MSS_dew(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
		double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF);
	int twoPhaseFlashEnv_BFGS_MSS2(double P, double P_old0, double T, double T_old, double* ai, double* bi, double** BIC, double* ZI,
		double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF);
private:
	int Nc;
	string EOS;
	/*double* Tc;
	double* Pc;
	double* omega;
	double* Mw;
	double* PARA;
	double* ZI;
	double** BIC;*/
};

