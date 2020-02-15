#pragma once



#include<string>
#include<iostream>


#include"global.h"
#include<cmath>



using namespace std;


class BlackOil
{
public:
	BlackOil();
	~BlackOil();
	BlackOil(double T, double P, double gammaGas, double API, double Rs, double Co,
		string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
		string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod);
	BlackOil(double gammaGas, double API, double Rs, double Co,
		string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
		string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod);
	BlackOil(double gammaGas, double API, double Rs, double Co, double waterSal,
		string units, string Pseudocritical_method, string MwMethod, string muGasMethod,
		string muDeadOilMethod, string muSatOilMethod, string ZFactorMethod, string BobMethod);
	double bubblePointPressure(string corr);
	double SaturatedGasOilRatio(string corr);
	double OFVF(string corr);
	void DeadOilViscosity(string corr);
	double SatOilViscosity();
	double GFVF();
	double GasViscosity(string corr);
	double PbStanding();
	double PbLasater();
	double PbGlaso();
	double PbBeggs_Vazquez();
	double PbDindoruk_Christman();
	double RsStanding();
	double RsLasater();
	double RsVazquez_Beggs();
	double RsGlaso();
	double RsDindoruk_Christman();
	double BobStanding();
	double BobStanding(double T, double Rs, double gammaGas, double gammaOil);
	double BobGlaso();
	double BobGlaso(double T, double Rs, double gammaGas, double gammaOil);
	double BobAlMarhoun();
	double BobAlMarhoun(double T, double Rs, double gammaGas, double gammaOil);
	double BobVazquez_Beggs();
	double BobVazquez_Beggs(double T, double Rs, double API, double gammaGas);
	double BobDindoruk_Christman();
	double BobDindoruk_Christman(double T, double Rs, double gammaGas, double gammaOil, double API);
	double BobFunc();
	double BoUnderSat();
	double muDeadBeal_Standing();
	double muDeadBeal_Standing(double T, double API);
	double muDeadBeggs_Robinson();
	double muDeadBeggs_Robinson(double T, double API);
	double muDeadGlaso();
	double muDeadGlaso(double T, double API);
	double muDeadAlKhafaji();
	double muDeadAlKhafaji(double T, double API);
	double muDeadDindoruk_Christman();
	double muDeadDindoruk_Christman(double T, double Pb, double API);
	double muDeadOilFunc();
	double muSatStanding();
	double muSatStanding(double Rs, double muDeadOil);
	double muSatBeggs_Robinson();
	double muSatBergman();
	double muSatAziz();
	double muSatAlKhafaji();
	double muSatDindoruk_Christman();
	double muSatOilFunc();
	double muUnderSatFunc();
	double BgCalculation();
	double BgCalculation(double P, double T, double Z);
	double ZFactor();
	void Pseudocritical_props(string corr);
	void Pc_Tc_Sutton();
	double muGasDempsey();
	double muGasLee_Gunzalez();
	double muGasLee_Gunzalez(double T, double MwGas, double dGas);
	double muGasLucas();
	double muGasLucas(double gasPpc, double gasTpc, double P, double T, double MwGas);
	double muGasFunc();
	double ZFactor_PVTsim();
	double ZFactor_PVTsim(double P, double T, double gasPpc, double gasTpc);
	double gammaOilCal();
	double gammaOilCal(double API);
	double stockTankMw(string corr);
	double GasDensity();
	double GasDensity(double P, double T, double Z);
	void gasCriticalProps(string corr);
	double gasMolecularWeight();
	double oilDensityBubble(double gammaOil, double Rs, double gammaGas, double Bo);
	double oilDensityBubble();
	double oilDensityAboveBubble();

	double muWaterFunc();
	double BwFunc();
	double waterDensityFunc();

	void set_T(double var);
	void set_P(double var);
	void set_Rs(double var);
	void set_Bob(double var);
	void set_Bg(double var);
	void set_gammaGas(double var);
	void set_gammaOil(double var);
	void set_API(double var);
	void set_Pb(double var);
	void set_MwOil(double var);
	void set_MwGas(double var);
	void set_Z(double var);
	void set_gasTpc(double var);
	void set_gasPpc(double var);
	void set_muDeadOil(double var);
	void set_muSatOil(double var);
	void set_muGas(double var);
	void set_muOil(double var);
	void set_dGas(double var);
	void set_dOil(double var);

	double get_T();
	double get_P();
	double get_Rs();
	double get_Bob();
	double get_Bg();
	double get_gammaGas();
	double get_gammaOil();
	double get_API();
	double get_Pb();
	double get_MwOil();
	double get_MwGas();
	double get_Z();
	double get_gasTpc();
	double get_gasPpc();
	double get_muDeadOil();
	double get_muSatOil();
	double get_muGas();
	double get_muOil();
	double get_dGas();
	double get_dOil();





private:

	double T;
	double P;
	double Rs;
	double Bob;
	double BoUnd;
	double Bg;
	double gammaGas;
	double gammaOil;
	double API;
	double Pb;
	double MwOil;
	double MwGas;
	double Z;
	double gasTpc;
	double gasPpc;
	double muDeadOil;
	double muSatOil;
	double muGas;
	double muOil;
	double dGas;
	double dOil;
	double dOilAboveBubble;
	double Co;
	double muUnsatOil;

	double muWater;
	double waterSal;
	double Bw;
	double rhoW;
	//// methods
	string muGasMethod;
	string muDeadOilMethod;
	string muSatOilMethod;
	string satOilVisMethod;
	string ZFactorMethod;
	string BobMethod;
};

