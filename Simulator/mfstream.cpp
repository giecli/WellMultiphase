
#include "mfstream.h"


mFstream::mFstream()
{
}


mFstream::~mFstream()
{
}


int mFstream::FileSearch(char *rSeek) {
	register int i;
	char str[MAX_STRING_LENGTH];

	clear();                 // clear fail and eof bits
	seekg(0, std::ios::beg); // back to the start!
	*str = '\0';
	do {
		i = ReadWord(str);
		if (!strcmp(str, rSeek)) return -1;
	} while (i);

	return 0;		//Nothing found
}

int mFstream::ReadWord(char *rWord) {
	char ch;
	register int i = 0;

	*rWord = '\0';
	do {
		get(ch);
		if (eof()) {
			return 0;		//Nothing has been read
		}
	} while ((ch < 33) || (ch > 126));

	while ((ch > 32) && (ch < 127)) {
		*(rWord + i) = ch;
		i++;
		get(ch);
		if (eof()) {
			*(rWord + i) = '\0';
			return 1;		//Read, but end of file also encountered
		}
	}

	*(rWord + i) = '\0';
	return -1;		//Correct execution
}

int mFstream::readPVT_DATAFILE(int& Nc, string address) {
	char str[MAX_STRING_LENGTH];
	mFstream pvtFile;
	pvtFile.open(address, ios::in);
	if (pvtFile.is_open()) {
		{char keyword[] = "NCOMPS";
		if (!pvtFile.FileSearch(keyword)) mException(F0010);
		if (!pvtFile.ReadWord(str)) mException(F0110);
		Nc = atoi(str); }
	}
	else {
		cout << F0001 << endl;
		cout << address << endl;
		return 1;
	}
	return 0;
}




int mFstream::readPVT_DATAFILE(int& Nc, double* Tc, double* Pc, string& EOS, vector<string>& compName,
	double* omega, double* Mw, double* PARA, double* ZI, double** BIC, double Tref, string address) {
	char str[MAX_STRING_LENGTH];
	mFstream pvtFile;
	pvtFile.open(address, ios::in);
	if (pvtFile.is_open()) {

		{char keyword[] = "EOS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		EOS = str; }

		{char keyword[] = "RTEMP";
		if (!pvtFile.FileSearch(keyword)) mException(F0030);
		if (!pvtFile.ReadWord(str)) mException(F0130);
		Tref = atof(str); }

		{char keyword[] = "CNAMES";
		if (!pvtFile.FileSearch(keyword)) mException(F0040);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0140);
			compName.push_back(str);
		}}

		{char keyword[] = "TCRIT";
		if (!pvtFile.FileSearch(keyword)) mException(F0050);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0150);
			Tc[i] = atof(str);
		}}

		{char keyword[] = "PCRIT";
		if (!pvtFile.FileSearch(keyword)) mException(F0060);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0160);
			Pc[i] = (1e+5) * atof(str);
		}}

		{char keyword[] = "ACF";
		if (!pvtFile.FileSearch(keyword)) mException(F0070);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0170);
			omega[i] = atof(str);
		}}

		{char keyword[] = "MW";
		if (!pvtFile.FileSearch(keyword)) mException(F0080);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0180);
			Mw[i] = 1e-3 * atof(str);
		}}

		{char keyword[] = "PARACHOR";
		if (!pvtFile.FileSearch(keyword)) mException(F0090);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0190);
			PARA[i] = atof(str);
		}}

		{char keyword[] = "ZI";
		if (!pvtFile.FileSearch(keyword)) mException(F0210);
		for (int i = 0; i < Nc; i++) {
			if (!pvtFile.ReadWord(str)) mException(F0310);
			ZI[i] = atof(str);
		}}

		{char keyword[] = "BIC";
		if (!pvtFile.FileSearch(keyword)) mException(F0220);
		for (int i = 0; i < Nc; i++) {
			for (int j = 0; j < i + 1; j++) {
				if (i == j) {
					BIC[i][j] = 0;
				}
				else {
					if (!pvtFile.ReadWord(str)) mException(F0320);
					BIC[i][j] = atof(str);
					BIC[j][i] = atof(str);
				}
			}
		}}
		return 0;
	}
	else {
		cout << F0001 << endl;
		cout << address << endl;
		return 1;
	}
}

int mFstream::readTPHFLUIDBLACKOILOILGAS_DATAFILE(double& WellDepth, double& theta, int& gridNum, double& d, 
	double& e, double& Qo, double& Qg,
	double& Pw, double& Rs, double& API, double& Co, double& Tw, double& dL, double& gammaGas, 
	string& Pseudocritical_method, string& MwMethod, string& muGasMethod,
	string& muDeadOilMethod, string& muSatOilMethod, string& ZFactorMethod, string& BobMethod, string& flowDir,
	string& waterProd, string& units, string& address){

	char str[MAX_STRING_LENGTH];
	mFstream pvtFile;
	pvtFile.open(address, ios::in);
	if (pvtFile.is_open()) {

		{char keyword[] = "WELLDEPTH";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		WellDepth = atof(str); }

		{char keyword[] = "INNERTUBINGD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		d = atof(str); }

		{char keyword[] = "THETA";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		theta = atof(str);
		theta = theta * piNum / 180.0; }

		{char keyword[] = "GRIDNUM";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		gridNum = atoi(str);
		}

		{char keyword[] = "ROUGHNESS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		e = atof(str); }

		{char keyword[] = "QOIL";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Qo = atof(str); }

		{char keyword[] = "PWELLHEAD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Pw = atof(str); }

		{char keyword[] = "TWELLHEAD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Tw = atof(str); }


		{char keyword[] = "SOLUTIONGAS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Rs = atof(str); }

		{char keyword[] = "API";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		API = atof(str); }

		{char keyword[] = "CO";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Co = atof(str); }

		{char keyword[] = "GAMMAGAS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		gammaGas = atof(str); }

		{char keyword[] = "Pseudocritical_method";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Pseudocritical_method = str; }

		{char keyword[] = "MwMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		MwMethod = str; }

		{char keyword[] = "muGasMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muGasMethod = str; }

		{char keyword[] = "muDeadOilMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muDeadOilMethod = str; }

		{char keyword[] = "muSatOilMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muSatOilMethod = str; }

		{char keyword[] = "ZFactorMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		ZFactorMethod = str; }

		{char keyword[] = "BobMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		BobMethod = str; }

		{char keyword[] = "flowDir";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		flowDir = str; }

		{char keyword[] = "UNIT";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		units = str; }

		if (units == "field") {
			WellDepth = foot_to_meter(WellDepth);
			e = inch_to_meter(e);
			d = inch_to_meter(d);
			Qo = STB_to_scm(Qo);
			Pw = psi_to_Pas(Pw);
			Tw = F_to_K(Tw);
			Co = Co / psi_to_Pas(1.0);
			Rs = scfPerSTB_to_scmPerscm(Rs);
		}

		dL = WellDepth / gridNum;
		Qg = Rs * Qo;
		return 0;
	}
	else {
		cout << F0001 << endl;
		cout << address << endl;
		return 1;
	}


	
}


int mFstream::readTPHFLUIDBLACKOILOILGASWATER_DATAFILE(double& WellDepth, double& theta, int& gridNum, double& d,
	double& e, double& Qo, double& Qg, double& Qw,
	double& Pw, double& Rs, double& API, double& Co, double& waterSal, double& Tw, double& dL, double& gammaGas,
	string& Pseudocritical_method, string& MwMethod, string& muGasMethod,
	string& muDeadOilMethod, string& muSatOilMethod, string& ZFactorMethod, string& BobMethod, string& flowDir,
	string& waterProd, string& units, string& address) {

	char str[MAX_STRING_LENGTH];
	mFstream pvtFile;
	pvtFile.open(address, ios::in);
	if (pvtFile.is_open()) {

		{char keyword[] = "WELLDEPTH";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		WellDepth = atof(str); }

		{char keyword[] = "INNERTUBINGD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		d = atof(str); }

		{char keyword[] = "THETA";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		theta = atof(str);
		theta = theta * piNum / 180.0; }

		{char keyword[] = "GRIDNUM";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		gridNum = atoi(str);
		}

		{char keyword[] = "ROUGHNESS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		e = atof(str); }

		{char keyword[] = "QOIL";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Qo = atof(str); }

		{char keyword[] = "QWATER";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Qw = atof(str); }

		{char keyword[] = "PWELLHEAD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Pw = atof(str); }

		{char keyword[] = "TWELLHEAD";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Tw = atof(str); }


		{char keyword[] = "SOLUTIONGAS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Rs = atof(str); }

		{char keyword[] = "API";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		API = atof(str); }

		{char keyword[] = "CO";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Co = atof(str); }

		{char keyword[] = "GAMMAGAS";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		gammaGas = atof(str); }

		{char keyword[] = "WATERSAL";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		waterSal = atof(str); }

		{char keyword[] = "Pseudocritical_method";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		Pseudocritical_method = str; }

		{char keyword[] = "MwMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		MwMethod = str; }

		{char keyword[] = "muGasMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muGasMethod = str; }

		{char keyword[] = "muDeadOilMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muDeadOilMethod = str; }

		{char keyword[] = "muSatOilMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		muSatOilMethod = str; }

		{char keyword[] = "ZFactorMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		ZFactorMethod = str; }

		{char keyword[] = "BobMethod";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		BobMethod = str; }

		{char keyword[] = "flowDir";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		flowDir = str; }

		{char keyword[] = "UNIT";
		if (!pvtFile.FileSearch(keyword)) mException(F0020);
		if (!pvtFile.ReadWord(str)) mException(F0120);
		units = str; }

		if (units == "field") {
			WellDepth = foot_to_meter(WellDepth);
			e = inch_to_meter(e);
			d = inch_to_meter(d);
			Qo = STB_to_scm(Qo);
			Qw = STB_to_scm(Qw);
			Pw = psi_to_Pas(Pw);
			Tw = F_to_K(Tw);
			Co = Co / psi_to_Pas(1.0);
			Rs = scfPerSTB_to_scmPerscm(Rs);
		}

		dL = WellDepth / gridNum;
		Qg = Rs * Qo;
		return 0;
	}
	else {
	cout << F0001 << endl;
	cout << address << endl;
	return 1;
	}
}

int mFstream::writeP_vs_L(int gridNum, double* P, double* L, string units, string caseName, string fileAddress) {
	mFstream writeObj;
	string fileName = "P_vs_L.txt";
	string saveAddress = fileAddress + "/" + fileName;
	writeObj.open(saveAddress, ios::out);
	if (writeObj.is_open()) {
		writeObj << "--Pressure versus Depth for " << caseName << endl;
		writeObj << "--Number of grids " << gridNum << endl;
		writeObj << "grid		Depth(ft)			Pressure(psi)" << endl;
		register int i;
		if ((units == "SI") || (units == "metric") || (units == "METRIC")) {
			writeObj << int(0 + 1) << "			" << L[0] << "					" << P[0] << endl;
			for (i = 1; i < gridNum; i++) {
				writeObj << int(i + 1) << "			" << L[i] << "			" << P[i] << endl;
			}
		}
		else if ((units == "Field") || (units == "field") || (units == "FIELD")) {
			writeObj << int(0 + 1) << "			" << m_to_foot(L[0]) << "				" << Pas_to_psi(P[0]) << endl;
			for (i = 1; i < gridNum; i++) {
				writeObj << int(i + 1) << "			" << m_to_foot(L[i]) << "			" << Pas_to_psi(P[i]) << endl;
			}
		}
		writeObj << "--End of file";
		writeObj.close();
		cout << "The results saved in :" << endl;
		cout << saveAddress << endl << endl;
		return 0;
	}
	else {
		cout << "Cannot open " << "P_vs_L.txt" << endl;
		return 1;
	}
}


int mFstream::writeTwoPhaseCompositionalFlashCalRes(int Nc, double* Zl, double* x, double* y, double T, double P,
	double nv, vector<string> components, string fileAddress) {

	mFstream writeObj;
	string fileName = "TwoPhaseFlashRes.txt";
	string saveAddress = fileAddress + "/" + fileName;
	writeObj.open(saveAddress, ios::out);
	if (writeObj.is_open()) {
		writeObj << "--Two phase flash calculation results " << endl;
		writeObj << "--Number of compoenents " << Nc << endl;
		writeObj << "Temperature (K)	" << T << endl;
		writeObj << "Pressure (Pas)	" << P << endl;

		writeObj << "Temperature (F)	" << K_to_F(T) << endl;
		writeObj << "Pressure (psi)	" << Pas_to_psi(P) << endl << endl;

		writeObj << "			    phase1		   phase2" << endl;
		writeObj << "mole%		   " << (1.0 - nv)*100 << "	      " << nv*100 << endl << endl;
		writeObj << "Nc    component     phase1            phase2" << endl;
		register int i;
		for (i = 0; i < Nc; i++) {
			writeObj << i << setw(10) << components[i] << setw(17) << x[i] << setw(17) << y[i] << endl;
		}
		writeObj << "--End of file";
		writeObj.close();
		cout << "The results saved in :" << endl;
		cout << saveAddress << endl << endl;
		return 0;
	}
	else {
		cout << "Cannot open " << "TwoPhaseFlashcalculationresults.txt" << endl;
		return 1;
	}
}