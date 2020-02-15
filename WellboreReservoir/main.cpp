#include"BlackOil.h"
#include"Compositional.h"
#include"twoPhaseModels.h"
#include"mFstream.h"
#include<iomanip>
#include<string>

#include <direct.h>
#define GetCurrentDir _getcwd






using namespace std;

void selectSimulation(int& simulation);
//void closeTheProgram(int & end);
string GetCurrentWorkingDir(void);



int main() {
	//twoPhaseModels pipe;
	mFstream readFileObj;
	//// Select pipe flow and compositional flash calculation
	//int simulation = -1;
	//int end = -1;
	//do {

	//	selectSimulation(simulation);
	//	closeTheProgram(end);

	//} while (end == 1);

	return 0;
}

void selectSimulation(int& simulation) {

	while ((simulation != 1) || (simulation != 2)) {
		
		cout << "Enter 1 to start flash calculation.\n" << endl;
		cout << "or \nEnter 0 to close the program." << endl;
		cin >> simulation;
		if (simulation == 1) {
			int Nc = 0;
			double* Tc;
			double* Pc;
			string EOS = "";
			vector<string> compName;
			double* omega;
			double* Mw;
			double* PARA;
			double* ZI;
			double** BIC;
			double Tref = 0.0;
			string fileAddress = GetCurrentWorkingDir();
			fileAddress.erase(fileAddress.end() - 17, fileAddress.end());
			fileAddress = fileAddress + "Examples";

			string fileName = "";
			cout << "Please enter file name." << endl;
			cin >> fileName;
			fileAddress = fileAddress + "\\" + fileName;

			//readFileObj.readPVT_DATAFILE(Nc, fileAddress);
			// initiallization
			Tc = new double[Nc];
			Pc = new double[Nc];
			omega = new double[Nc];
			Mw = new double[Nc];
			PARA = new double[Nc];
			ZI = new double[Nc];
			BIC = new double*[Nc];
			for (int i = 0; i < Nc; i++) {
				BIC[i] = new double[Nc];
			}
			//readFileObj.readPVT_DATAFILE(Nc, Tc, Pc, EOS, compName, omega, Mw, PARA, ZI, BIC, Tref, fileAddress);
		}
		else if (simulation == 0) {
			cout << "quit? Y or N" << endl;
			string endCommand = "";
			cin >> endCommand;
			if (endCommand == "Y" || endCommand == "y") {
				break;
			}
		}
		else {
			cout << "Wrong selection./nPlease note: " << endl;
		}
	}
}
//
//void closeTheProgram(int & end) {
//	cout << "Simulation ended successfully." << endl;
//	while ((end != 0) || (end != 1))
//	{
//		cout << "Enter 0 to close the program or Enter 1 to start a new simulation." << endl;
//		cin >> end;
//		if ((end != 0) || (end != 1)) {
//			cout << "Wrong selection.\nPlease note: " << endl;
//			break;
//		}
//	}
//}
//
//
string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	string current_working_dir(buff);
	return current_working_dir;
}