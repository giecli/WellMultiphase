#pragma once
#include<fstream>
#include"mException.h"
#include"ErrorCode.h"
#include<vector>
#include"global.h"
#include<string>
#include<iomanip>

#define MAX_STRING_LENGTH 256
using namespace std;

class mFstream : public fstream
{
public:
	mFstream();
	~mFstream();
	int FileSearch(char*);
	int ReadWord(char*);

	//// reading compositional data
	int readPVT_DATAFILE(int& Nc, string address);
	int readPVT_DATAFILE(int& Nc, double* Tc, double* Pc, string& EOS, vector<string>& compName,
		double* omega, double* Mw, double* PARA, double* ZI, double** BIC, double Tref, string address);
	/////////////////////////////////***************************////////////////////////////////////////

	/////////
	int readTPHFLUIDBLACKOILOILGAS_DATAFILE(double& WellDepth, double& theta, int& gridNum, 
		double& d, double& e, double& Qo, double& Qg, double& Pw, double& Rs,
		double& API, double& Co, double& Tw, double& dL, double& gammaGas,
		string& Pseudocritical_method, string& MwMethod, string& muGasMethod,
		string& muDeadOilMethod, string& muSatOilMethod, string& ZFactorMethod, string& BobMethod,
		string& flowDir, string& waterProd, string& units, string& address);
	int readTPHFLUIDBLACKOILOILGASWATER_DATAFILE(double& WellDepth, double& theta, int& gridNum,
		double& d, double& e, double& Qo, double& Qg, double& Qw, double& Pw, double& Rs,
		double& API, double& Co, double& waterSal, double& Tw, double& dL, double& gammaGas,
		string& Pseudocritical_method, string& MwMethod, string& muGasMethod,
		string& muDeadOilMethod, string& muSatOilMethod, string& ZFactorMethod, string& BobMethod,
		string& flowDir, string& waterProd, string& units, string& address);

	//// write pressure vs depth to text
	int writeP_vs_L(int gridNum, double* P, double* L, string unit, string caseName, string fileAddress);
	
	//// write compositional flash calculation to text
	int writeTwoPhaseCompositionalFlashCalRes(int Nc, double* Zl, double* x, double* y, double T, double P,
		double nv, vector<string> components, string fileAddress);

};

