
#include"mfstream.h"
#include"BlackOil.h"
#include"Compositional.h"
#include"twoPhaseModels.h"

#include<iomanip>
#include<string>

#include <direct.h>
#define GetCurrentDir _getcwd

#include <stdio.h>  // defines FILENAME_MAX







using namespace std;

void selectSimulation(int& simulation);
void closeTheProgram(int & end);
string GetCurrentWorkingDir(void);
string GetCurrentWorkingDir1();
void twoPhaseFlashSimulation();
void twoPhaseFluidFlowInPipeSimulation_OILGAS();
void twoPhaseFluidFlowInPipeSimulation_OILGASWATER();

int main() {
	twoPhaseModels pipe;
	
	//// Select pipe flow and compositional flash calculation
	int simulation = -1;
	int end = -1;
	do {

		selectSimulation(simulation);
		closeTheProgram(end);

	} while (end == 1);

	return 0;
}

void selectSimulation(int& simulation) {

	while ((simulation != 1) || (simulation != 2) || (simulation != 3)) {

		cout << "Enter 1 to start flash calculation." << endl;
		cout << "Enter 2 to start steady state two phase (oil and gas) flow in the well." << endl;
		cout << "Enter 3 to start steady state two phase (oil, gas and water) flow in the well." << endl;
		cout << "or \nEnter 0 to close the program." << endl;
		cin >> simulation;
		if (simulation == 1) {
			twoPhaseFlashSimulation();
			// write results in the text
		}
		else if (simulation == 0) {
			cout << "quit? Y or N" << endl;
			string endCommand = "";
			cin >> endCommand;
			if (endCommand == "Y" || endCommand == "y") {
				break;
			}
		}
		else if (simulation == 2) {
			twoPhaseFluidFlowInPipeSimulation_OILGAS();
			// write results in the text
		}
		else if (simulation == 3) {
			twoPhaseFluidFlowInPipeSimulation_OILGASWATER();
			// write results in the text
		}
		else {
			cout << "Wrong selection./nPlease note: " << endl;
		}
	}
}

void closeTheProgram(int & end) {
	cout << "Simulation ended successfully." << endl;
	while ((end != 0) || (end != 1))
	{
		cout << "Enter 0 to close the program or Enter 1 to start a new simulation." << endl;
		cin >> end;
		if ((end != 0) || (end != 1)) {
			cout << "Wrong selection.\nPlease note: " << endl;
			break;
		}
	}
}


string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	string current_working_dir(buff);
	return current_working_dir;
}

std::string GetCurrentWorkingDir1()
{
	std::string cwd("\0", FILENAME_MAX );
	return _getcwd(&cwd[0], cwd.capacity());
}

void twoPhaseFlashSimulation() {
	int Nc = 0;
	double* Tc;
	double* Pc;
	string EOS = "";
	vector<string> compName;
	double* ACF;
	double* Mw;
	double* PARA;
	double* ZI;
	double** BIC;
	double Tref = 0.0;
	string fileAddress = GetCurrentWorkingDir1();
	fileAddress.erase(fileAddress.end() - 9, fileAddress.end());
	fileAddress = fileAddress + "Examples";

	string fileName = "";
	cout << "Please enter file name." << endl;
	cin >> fileName;
	fileAddress = fileAddress + "\\" + fileName;
	mFstream readFileObj;
	string saveAddress;
	cout << "Please enter a directory to save the results." << endl;
	cout << "For example: D:/PVTMODEL/well1" << endl;
	cin >> saveAddress;

	int suc = readFileObj.readPVT_DATAFILE(Nc, fileAddress);
	// initiallization
	if (suc == 0) {
		Tc = new double[Nc];
		Pc = new double[Nc];
		ACF = new double[Nc];
		Mw = new double[Nc];
		PARA = new double[Nc];
		ZI = new double[Nc];
		BIC = new double*[Nc];
		for (int i = 0; i < Nc; i++) {
			BIC[i] = new double[Nc];
		}
		double* ai = new double[Nc];
		double* bi = new double[Nc];

		readFileObj.readPVT_DATAFILE(Nc, Tc, Pc, EOS, compName, ACF, Mw, PARA, ZI, BIC, Tref, fileAddress);
		if (Nc > 1) {
			BIC[0][1] = -BIC[0][1];
			BIC[1][0] = -BIC[1][0];
		}
		double P, T;
		cout << "Enter Pressure (bar):" << endl;
		cin >> P;
		P = bar_to_Pas(P);
		cout << "Enter Temperature (C):" << endl;
		cin >> T;
		T = C_to_K(T);
		Compositional compObj(Tc, Pc, ACF, ai, bi, Nc, T, P, EOS);
		double* K = new double[Nc];
		double* x = new double[Nc];
		double* y = new double[Nc];
		double nv = 0.5;
		compObj.wilsonK(P, T, Tc, Pc, ACF, K);
		double sumU = compObj.twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, ZI, K, nv);
		//compObj.twoPhaseStabilityFlash(P, T, ai, bi, BIC, ZI, K, nv);
		if (sumU > 1.0) {
			cout << "Enter first guess for nv:" << endl;
			cin >> nv;
			compObj.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, ZI, K, nv, x, y);
			cout << "nv is: " << nv << endl;

			mFstream writeObj;
			writeObj.writeTwoPhaseCompositionalFlashCalRes(Nc, ZI, x, y, T, P, nv, compName, saveAddress);
		}
		else {
			cout << "The fluid at " << Pas_to_bar(P) << " (bar)	and " << K_to_C(T) << " (C)	is not two phase." << endl;
		}
		delete[] Tc;
		delete[] Pc;
		delete[] ACF;
		delete[] Mw;
		delete[] PARA;
		delete[] ZI;
		register int i;
		for (i = 0; i < Nc; i++) {
			delete[] BIC[i];
		}
		delete[] BIC;
		delete[] ai;
		delete[] bi;

		delete[] K;
		delete[] x;
		delete[] y;
	}
	else if (suc == 1) {
		string errmsg = fileAddress + " : doesn't exist.";
		cout << errmsg << endl << endl;
	}
	
	
	
}

void twoPhaseFluidFlowInPipeSimulation_OILGAS() {

	string fileAddress = GetCurrentWorkingDir1();
	fileAddress.erase(fileAddress.end() - 9, fileAddress.end());
	fileAddress = fileAddress + "Examples";

	string fileName = "";
	cout << "Please enter file name." << endl;
	cin >> fileName;
	fileAddress = fileAddress + "\\" + fileName;
	mFstream readFileObj;
	string saveAddress;
	cout << "Please enter a directory to save the results." << endl;
	cout << "For example: D:/wellModels/well1" << endl;
	cin >> saveAddress;

	double WellDepth; 	double theta; 	double d;	double e;
	double Qo;	double Pw; double Rs; double API; double Co;
	double Tw; double dL; double gammaGas; int gridNum; double Qg;
	string Pseudocritical_method; string MwMethod; string muGasMethod;
	string muDeadOilMethod; string muSatOilMethod; string ZFactorMethod;
	string BobMethod; string flowDir; string waterProd; string units;
	
	double* PAns; double* L;
		
	int suc = readFileObj.readTPHFLUIDBLACKOILOILGAS_DATAFILE(WellDepth, theta, gridNum, d, e,
		Qo, Qg, Pw, Rs, API, Co, Tw, dL, gammaGas, Pseudocritical_method, MwMethod, muGasMethod,
		muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir,
		waterProd, units, fileAddress);
	if (suc == 0) {
		twoPhaseModels FlowProfileObj;
		PAns = new double[gridNum + 1];
		L = new double[gridNum + 1];
		FlowProfileObj.FlowProfile(Qg, Qo, Pw, gridNum, Tw, Pw, gammaGas, API,
			Rs, Co, e, d, dL, theta, "SI", Pseudocritical_method,
			MwMethod, muGasMethod, muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir,
			PAns, L);

		mFstream writeObj;
		string caseName = "Two phase flow (oil and gas) in the well";
		writeObj.writeP_vs_L(gridNum + 1, PAns, L, units, caseName, saveAddress);
		delete[] PAns;
		delete[] L;
	}
	else if (suc == 1) {
		string errmsg = fileAddress + " : doesn't exist.";
		cout << errmsg << endl << endl;
	}

}


void twoPhaseFluidFlowInPipeSimulation_OILGASWATER() {

	string fileAddress = GetCurrentWorkingDir1();
	fileAddress.erase(fileAddress.end() - 9, fileAddress.end());
	fileAddress = fileAddress + "Examples";

	string fileName = "";
	cout << "Please enter file name." << endl;
	cin >> fileName;
	fileAddress = fileAddress + "\\" + fileName;
	mFstream readFileObj;
	string saveAddress;
	cout << "Please enter a directory to save the results." << endl;
	cout << "For example: D:/wellModels/well1" << endl;
	cin >> saveAddress;

	double WellDepth; 	double theta; 	double d;	double e;
	double Qo;	double Pw; double Rs; double API; double Co;
	double Tw; double dL; double gammaGas; int gridNum; double Qg;
	double Qw; double waterSal;
	string Pseudocritical_method; string MwMethod; string muGasMethod;
	string muDeadOilMethod; string muSatOilMethod; string ZFactorMethod;
	string BobMethod; string flowDir; string waterProd; string units;

	double* PAns; double* L;

	int suc = readFileObj.readTPHFLUIDBLACKOILOILGASWATER_DATAFILE(WellDepth, theta, gridNum, d, e,
		Qo, Qg, Qw, Pw, Rs, API, Co, waterSal, Tw, dL, gammaGas, Pseudocritical_method, MwMethod, muGasMethod,
		muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir,
		waterProd, units, fileAddress);
	if (suc == 0) {
		PAns = new double[gridNum + 1];
		L = new double[gridNum + 1];
		twoPhaseModels FlowProfileObj;
		FlowProfileObj.FlowProfile(Qg, Qo, Qw, waterSal, Pw, gridNum, Tw, Pw, gammaGas, API,
			Rs, Co, e, d, dL, theta, "SI", Pseudocritical_method,
			MwMethod, muGasMethod, muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir, PAns, L);
		mFstream writeObj;
		string caseName = "Two phase flow (oil, gas and water) in the well";
		writeObj.writeP_vs_L(gridNum + 1, PAns, L, units, caseName, saveAddress);
		delete[] PAns;
		delete[] L;
	}
	else if (suc == 1) {
		string errmsg = fileAddress + " : doesn't exist.";
		cout << errmsg << endl;
	}
	
}