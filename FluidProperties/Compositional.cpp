#include "Compositional.h"

Compositional::Compositional(double* Tc, double* Pc, double* ACF, double* ai, double* bi,
	int inNc, double T, double P, string inEOS)
{
	Nc = inNc;
	EOS = inEOS;

	if (EOS == "RK") {
		RKEOS(T, P, Tc, Pc, ai, bi);
	}
	else if (EOS == "SRK") {
		SRKEOS(T, P, Tc, Pc, ACF, ai, bi);
	}
	if (EOS == "GD") {
		GDEOS(T, P, Tc, Pc, ACF, ai, bi);
	}
	else if (EOS == "PR") {
		PREOS(T, P, Tc, Pc, ACF, ai, bi);
	}
	else if ("PR78") {
		PR78EOS(T, P, Tc, Pc, ACF, ai, bi);
	}


	/*double* Tc = new double[Nc];
	double* Pc = new double[Nc];
	double* omega = new double[Nc];
	double* Mw = new double[Nc];
	double* PARA = new double[Nc];
	double* ZI = new double[Nc];
	double** BIC = new double*[Nc];
	for (int i = 0; i < Nc; i++) {
		BIC[i] = new double[Nc];
	}

	for (int i = 0; i < Nc; i++) {
		Tc[i] = inTc[i];
		Pc[i] = inPc[i];
		omega[i] = inomega[i];
		Mw[i] = in
	}
*/
}


Compositional::~Compositional()
{
}
Compositional::Compositional()
{
}

double Compositional::cubicEOS(double T, double P, double* comp, double* ai, double* bi,
	double& a, double& b, double* S, double& A, double& B, double** BIC, double& Z, string& phase) {

	register int i, j;
	a = 0;
	A = 0;
	b = 0;
	B = 0;
	double c;
	for (i = 0; i < Nc; i++) {
		b += comp[i] * bi[i];
		for (j = 0; j < Nc; j++) {
			S[i] += comp[j] * (1 - BIC[i][j]) * sqrt(ai[j]);
		}
		S[i] *= sqrt(ai[i]);
		a += comp[i] * S[i];
	}
	if ((EOS == "RK") || (EOS == "SRK") || (EOS == "GD")) {
		c = 0;
	}
	else if ((EOS == "PR") || ("PR78")) {
		c = 1.0;
	}
	A = a * P / (Rgas * Rgas * T * T);
	B = b * P / (Rgas * T);


	mathFunc obj;
	double a1 = -(1 - c * B);
	double a2 = A - B * (1 + c) - B * B * (1 + 2 * c);
	double a3 = -(A * B - c * (B * B * B + B * B));
	vector<double> ans0, ans;
	obj.cubicEqSolver(a1, a2, a3, ans0);
	for (int i = 0; i < ans0.size(); i++) {
		if (ans0[i] > 0) {
			ans.push_back(ans0[i]);
		}
	}

	if (ans.size() == 1) {
		Z = ans[0];
		return Z;
	}
	else if (ans.size() == 2) {
		double ZA = ans[0];
		double ZB = ans[1];
		Z = cubicRootSelection(ZA, ZB, B, A, phase);
		return Z;
	}
	else if (ans.size() == 3) {
		double ZA, ZB;
		if (((ans[0] > ans[1]) && (ans[0] < ans[2])) || ((ans[0] < ans[1]) && (ans[0] > ans[2]))) {
			ZA = ans[1];
			ZB = ans[2];
		}
		else if (((ans[1] > ans[0]) && (ans[1] < ans[2])) || ((ans[1] < ans[0]) && (ans[1] > ans[2]))) {
			ZA = ans[0];
			ZB = ans[2];
		}
		else if (((ans[2] > ans[0]) && (ans[2] < ans[1])) || ((ans[2] < ans[0]) && (ans[2] > ans[1]))) {
			ZA = ans[0];
			ZB = ans[1];
		}
		Z = cubicRootSelection(ZA, ZB, B, A, phase);
		return Z;
	}
	/////// condition????
}

void Compositional::RKEOS(double T, double P, double* Tc, double* Pc, double* a, double* b) {
	double ohmA, ohmB;
	ohmA = 0.42747;
	ohmB = 0.08664;
	register int i;
	for (i = 0; i < Nc; i++) {
		a[i] = (ohmA / sqrt(T / Tc[i]) * Rgas * Rgas * Tc[i] * Tc[i] / Pc[i]);
		b[i] = (ohmB * Rgas * Tc[i] / Pc[i]);
	}

}
void Compositional::SRKEOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b) {
	double ohmA, ohmB;
	ohmA = 0.42747;
	ohmB = 0.08664;
	double alpha, m;
	register int i;
	for (i = 0; i < Nc; i++) {
		m = 0.480 + 1.574 * ACF[i] - 0.176 * ACF[i] * ACF[i];
		alpha = (1 + m * (1 - sqrt(T / Tc[i]))) * (1 + m * (1 - sqrt(T / Tc[i])));
		a[i] = (ohmA * alpha * Rgas * Rgas * Tc[i] * Tc[i] / Pc[i]);
		b[i] = (ohmB * Rgas * Tc[i] / Pc[i]);
	}
}
void Compositional::GDEOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b) {
	double ohmA, ohmB;
	ohmA = 0.42747;
	ohmB = 0.08664;
	double alpha, m;
	register int i;
	for (i = 0; i < Nc; i++) {
		m = 0.48508 + +1.55171 * ACF[i] - 0.15613 * ACF[i] * ACF[i];
		alpha = (1 + m * (1 - sqrt(T / Tc[i]))) * (1 + m * (1 - sqrt(T / Tc[i])));
		a[i] = (ohmA * alpha * Rgas * Rgas * Tc[i] * Tc[i] / Pc[i]);
		b[i] = (ohmB * Rgas * Tc[i] / Pc[i]);
	}
}
void Compositional::PREOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b) {
	double ohmA, ohmB;
	ohmA = 0.457235;
	ohmB = 0.077796;
	double alpha, m;
	register int i;
	for (i = 0; i < Nc; i++) {
		m = 0.37464 + 1.5422 * ACF[i] - 0.269926 * ACF[i] * ACF[i];
		alpha = (1 + m * (1 - sqrt(T / Tc[i]))) * (1 + m * (1 - sqrt(T / Tc[i])));
		a[i] = (ohmA * alpha * Rgas * Rgas * Tc[i] * Tc[i] / Pc[i]);
		b[i] = (ohmB * Rgas * Tc[i] / Pc[i]);
	}
}
void Compositional::PR78EOS(double T, double P, double* Tc, double* Pc, double* ACF, double* a, double* b) {
	double ohmA, ohmB;
	ohmA = 0.457235;
	ohmB = 0.077796;
	double alpha, m;
	register int i;
	for (i = 0; i < Nc; i++) {
		if (ACF[i] <= 0.49) {
			m = 0.37464 + 1.5422 * ACF[i] - 0.269926 * ACF[i] * ACF[i];
		}
		else
		{
			m = 0.379642 + ACF[i] * (1.48503 - 0.164423 * ACF[i] + +0.01666 * ACF[i] * ACF[i]);
		}
		alpha = (1 + m * (1 - sqrt(T / Tc[i]))) * (1 + m * (1 - sqrt(T / Tc[i])));
		a[i] = (ohmA * alpha * Rgas * Rgas * Tc[i] * Tc[i] / Pc[i]);
		b[i] = (ohmB * Rgas * Tc[i] / Pc[i]);
	}
}


void Compositional::cubicFugacityCoef(double A, double B, double a, double b, double* bi,
	double* ai, double* Si, double Z, double* phi) {
	register int i;
	double temp1, temp2, temp3;
	double d1, d2;
	if ((EOS == "RK") || (EOS == "SRK") || (EOS == "GD")) {
		d1 = 0.0;
		d2 = 1.0;
	}
	else if ((EOS == "PR") || ("PR78")) {
		d1 = 1 - sqrt(2);
		d2 = 1 + sqrt(2);
	}
	for (i = 0; i < Nc; i++) {
		temp1 = (bi[i] / b) * (Z - 1.0) - log(Z - B);
		temp2 = -1.0 / (d2 - d1) * (A / B) * (2 * Si[i] / a - bi[i] / b);
		temp3 = log((Z + d2 * B) / (Z + d1 * B));
		phi[i] = exp(temp1 + temp2 * temp3);
	}
}

void Compositional::wilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* K) {
	register int i;
	for (i = 0; i < Nc; i++) {
		K[i] = (Pc[i] / P) * exp(5.37 * (1 + ACF[i]) * (1 - Tc[i] / T));
	}
}


void Compositional::dTwilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* dK) {

	register int i;
	for (i = 0; i < Nc; i++) {
		dK[i] = (Pc[i] / P) * exp(5.37 * (1 + ACF[i]) * (1 - Tc[i] / T)) * ((Tc[i] / T / T) * 5.37 * (1 + ACF[i]));
	}
}


void Compositional::dPwilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* dK) {

	register int i;
	for (i = 0; i < Nc; i++) {
		dK[i] = (Pc[i] / P) * exp(5.37 * (1 + ACF[i]) * (1 - Tc[i] / T)) * (-1.0 / P);
	}
}


//void Compositional::mWilsonK(double P, double T, double* Tc, double* Pc, double* ACF, double* K) {
//	register int i;
//	for (i = 0; i < Nc; i++) {
//		K[i] = (Pc[i] / P) * exp(5.37 * (1 + ACF[i]) * (1 - Tc[i] / T));
//	}
//}

double Compositional::cubicRootSelection(double ZA, double ZB, double B, double A, string& phase) {
	double d1, d2;
	if ((EOS == "RK") || (EOS == "SRK") || (EOS == "GD")) {
		d1 = 0.0;
		d2 = 1.0;
	}
	else if ((EOS == "PR") || ("PR78")) {
		d1 = 1 - sqrt(2);
		d2 = 1 + sqrt(2);
	}
	double temp1 = log((ZB - B) / (ZA - B)) - (ZB - ZA);
	double temp2 = (1.0 / (d2 - d1)) * (A / B) * log(((ZB + d2 * B) * (ZA + d1 * B)) / (ZA + d2 * B) * (ZB + d1 * B));
	double term = temp1 + temp2;

	if (term > 0) {
		if (ZB > ZA) {
			phase = "VAPOR";
		}
		else {
			phase = "LIQUID";
		}
		return ZB;
	}
	else {
		if (ZA > ZB) {
			phase = "VAPOR";
		}
		else {
			phase = "LIQUID";
		}
		return ZA;
	}

}

int Compositional::twoPhaseFlash_SS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	register int i;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];

	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		K_old[i] = 0.0;
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}

	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	double f, df;
	int outer_counter = 0;
	int inner_counter = 0;
	do {
		outer_counter++;
		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);
		do {

			f = 0.0;
			df = 0.0;
			inner_counter = 0;
			for (i = 0; i < Nc; i++) {
				inner_counter++;
				f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
				df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
			}
			nv_old = nv;
			nv -= f / df;
			err_2 = abs(nv - nv_old);
			if (nv > 1.0 / (1.0 - Kmin)) {
				nv = 1.0 / (1.0 - Kmin);
			}
			if (nv < 1.0 / (1.0 - Kmax)) {
				nv = 1.0 / (1.0 - Kmax);
			}
		} while (err_2 > tol_2);

		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			K_old[i] = K[i];
			fL[i] = x[i] * P * phiL[i];
			fG[i] = y[i] * P * phiG[i];
			K[i] = K_old[i] * (fL[i] / fG[i]); // K_old[i] * (fL[i] / fG[i]); // phiL[i] / phiG[i];
		}
		err_1 = 0;
		for (i = 1; i < Nc; i++) {
			err_1 += (1 - fL[i] / fG[i]) * (1 - fL[i] / fG[i]);
		}

	} while (err_1 > tol_1);


	delete[] K_old;
	return 0;
}

int Compositional::twoPhaseFlash_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	register int i;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];

	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		K_old[i] = 0.0;
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}

	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	double f, df;
	int outer_counter = 0;
	int inner_counter = 0;
	do {
		outer_counter++;
		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= f / df;
		err_2 = abs(nv - nv_old);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 1.0 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.0 / (1.0 - Kmax);
		}

		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			K_old[i] = K[i];
			fL[i] = x[i] * P * phiL[i];
			fG[i] = y[i] * P * phiG[i];
			K[i] = K_old[i] * (fL[i] / fG[i]); // K_old[i] * (fL[i] / fG[i]); // phiL[i] / phiG[i];
		}
		err_1 = 0;
		for (i = 1; i < Nc; i++) {
			err_1 += (1 - fL[i] / fG[i]) * (1 - fL[i] / fG[i]);
		}

	} while (err_1 > tol_1);


	delete[] K_old;
	return 0;
}

/*
int Compositional::twoPhaseFlash_QNN_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	mathFunc mathSolver;
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_2 = 1.0;
	double nv_old;
	double* LnK_ = new double[Nc];
	double* LnK_0 = new double[Nc];
	double* LnK_1 = new double[Nc];
	double* g0 = new double[Nc];
	double* g1 = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double* temp1 = new double[Nc];
	double* temp2 = new double[Nc];
	double* temp3 = new double[Nc];
	double* y1 = new double[Nc];
	double* s1 = new double[Nc];
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 0; i < Nc; i++) {
		LnK_1[i] = log(K[i]);
		LnK_0[i] = log(K[i]);
		LnK_[i] = log(K[i]);
		g0[i] = 0.0;
		g1[i] = 0.0;
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
		temp1[i] = 0.0;
		temp2[i] = 0.0;
		temp3[i] = 0.0;
		y1[i] = 0.0;
		s1[i] = 0.0;
	}


	double f, df;
	int outer_counter = 0;
	int inner_counter = 0;
	double** H0 = new double*[Nc];
	double** H1 = new double*[Nc];
	double** mTemp1 = new double*[Nc];
	double** mTemp2 = new double*[Nc];
	double sTemp = 0.0;
	for (i = 0; i < Nc; i++) {
		H0[i] = new double[Nc];
		H1[i] = new double[Nc];
		mTemp1[i] = new double[Nc];
		mTemp2[i] = new double[Nc];
	}

	double err_1;



	do {
		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				H0[i][j] = 0.0;
				H1[i][j] = 0.0;
				mTemp1[i][j] = 0.0;
				mTemp2[i][j] = 0.0;
			}
		}

		for (i = 0; i < Nc; i++) {
			H0[i][i] = 1.0;
		}
		mathSolver.assignVecToVec(LnK_0, LnK_, Nc);
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);
		for (i = 0; i < Nc; i++) {
			g0[i] = LnK_0[i] - log(phiL[i]) + log(phiG[i]);
		}
		//mathSolver.assignVecToVec(g1, g0, Nc);
		mathSolver.matrixA_product_vectorb(H0, g0, temp1, Nc);
		mathSolver.vectorSubtraction(LnK_0, temp1, LnK_1, Nc);
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK_1[i]);
		}
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		do {
			//mathSolver.assignVecToVec(LnK_1, LnK_0, Nc);
			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);
			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);
			for (i = 0; i < Nc; i++) {
				g1[i] = LnK_1[i] - log(phiL[i]) + log(phiG[i]);
			}

			mathSolver.vectorSubtraction(g1, g0, y1, Nc);
			mathSolver.vectorSubtraction(LnK_1, LnK_0, s1, Nc);
			mathSolver.matrixA_product_vectorb(H0, g1, temp1, Nc);
			mathSolver.vectorb_product_matrixA(H0, s1, temp2, Nc);
			mathSolver.vectorA_CrossProduct_VectorAT(temp1, temp2, mTemp1, Nc);
			sTemp = mathSolver.vectorDotProduct(temp2, y1, Nc);
			mathSolver.matrixDividedByNumber(mTemp1, sTemp, Nc);
			mathSolver.matrixSummation(H0, mTemp1, H1, Nc);
			mathSolver.assignMatToMat(H1, H0, Nc);
			mathSolver.assignVecToVec(g1, g0, Nc);
			mathSolver.matrixA_product_vectorb(H0, g0, temp1, Nc);
			mathSolver.assignVecToVec(LnK_1, LnK_0, Nc);
			mathSolver.vectorSubtraction(LnK_1, temp1, LnK_1, Nc);
			err_1 = mathSolver.errorSolution(LnK_1, LnK_0, Nc);

			for (i = 0; i < Nc; i++) {
				K[i] = exp(LnK_1[i]);
			}
			for (i = 0; i < Nc; i++) {
				x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
				y[i] = K[i] * x[i];
			}
			normalizArr(x, Nc);
			normalizArr(y, Nc);
		} while (err_1 > 1e-5);
		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= f / df;

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);
		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 1.0 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.0 / (1.0 - Kmax);
		}

		outer_counter++;
		err_1 = mathSolver.errorSolution(LnK_1, LnK_, Nc);
	} while (err_1 > 1e-5);

	return 0;
}
*/

int Compositional::twoPhaseFlash_QNN_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-18;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc), g0(Nc), LnK1(Nc), Lnk0(Nc), LnphiL(Nc),
		s1(Nc), s1t(Nc), y1(Nc), LnphiG(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);

	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnK1(i) = log(K[i]);
	}

	g1 = LnK1 - LnphiL + LnphiG;

	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		Lnk0 = LnK1;

		LnK1 = Lnk0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);
		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}

		g1 = LnK1 - LnphiL + LnphiG;

		s1 = LnK1 - Lnk0;

		y1 = g1 - g0;
		H1 = H0 + (s1 - H0 * y1) * (s1.transpose() * H0);
		MatrixXd temp = (s1.transpose() * H0 * y1);

		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				H1(i, j) = H1(i, j) / temp(0, 0);
			}
		}

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= f / df;
		err_2 = abs(nv - nv_old);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 1.0 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.0 / (1.0 - Kmax);
		}

		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
	}

	return 0;
}


int Compositional::twoPhaseFlash_BFGS_MSS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc), g0(Nc), LnK1(Nc), Lnk0(Nc), LnphiL(Nc),
		s1(Nc), s1t(Nc), y1(Nc), LnphiG(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);

	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnK1(i) = log(K[i]);
	}

	g1 = LnK1 - LnphiL + LnphiG;

	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		Lnk0 = LnK1;

		LnK1 = Lnk0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);
		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}

		g1 = LnK1 - LnphiL + LnphiG;

		s1 = LnK1 - Lnk0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= f / df;
		err_2 = abs(nv - nv_old);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 1.0 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.0 / (1.0 - Kmax);
		}

		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
	}

	return 0;
}

int Compositional::twoPhaseFlash_BFGS_FullyImplicit(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv, double* x, double* y) {
	register int i, j, ii;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc + 1), g0(Nc + 1), LnK1(Nc + 1), Lnk0(Nc + 1), LnphiL(Nc),
		s1(Nc + 1), s1t(Nc + 1), y1(Nc + 1), LnphiG(Nc);
	MatrixXd H1(Nc + 1, Nc + 1), H0(Nc + 1, Nc + 1);

	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc + 1, Nc + 1);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnK1(i) = log(K[i]);
	}
	for (i = 0; i < Nc; i++) {
		g1(i) = LnK1(i) + LnphiL(i) - LnphiG(i);
	}
	f = 0.0;
	for (i = 0; i < Nc; i++) {
		f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
	}
	LnK1(Nc) = nv;
	g1(Nc) = f;

	double dphiLdx, dphiGdy, dLnphiLdLnK, dLnphiGdLnK;
	double* tempx = new double[Nc];
	double* tempy = new double[Nc];
	double* dphiLdx_f = new double[Nc];
	double* dphiLdx_b = new double[Nc];
	double* dphiGdy_f = new double[Nc];
	double* dphiGdy_b = new double[Nc];
	for (ii = 0; ii < Nc; ii++) {
		tempx[ii] = x[ii];
		tempy[ii] = y[ii];
	}
	for (i = 0; i < Nc; i++) {
		for (j = 0; j < Nc; j++) {

			tempx[j] = x[j] * 1.001;
			tempy[j] = y[j] * 1.001;
			normalizArr(tempx, Nc);
			normalizArr(tempy, Nc);
			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, tempx, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, dphiLdx_f);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, dphiGdy_f);

			for (ii = 0; ii < Nc; ii++) {
				tempx[ii] = x[ii];
				tempy[ii] = y[ii];
			}

			tempx[j] = x[j] * 0.999;
			tempy[j] = y[j] * 0.999;
			normalizArr(tempx, Nc);
			normalizArr(tempy, Nc);
			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, tempx, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, dphiLdx_b);
			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, dphiGdy_b);
			dphiLdx = (dphiLdx_f[j] - dphiLdx_b[j]) / (2 * 0.001 * x[j]);
			dphiGdy = (dphiGdy_f[j] - dphiGdy_b[j]) / (2 * 0.001 * y[j]);
			dLnphiLdLnK = -y[j] / K[i] / LnphiL[i] * dphiLdx;
			dLnphiGdLnK = x[j] * K[i] / LnphiG[i] * dphiGdy;

			if (i == j) {
				H1(i, j) = 1 + dLnphiLdLnK - dLnphiGdLnK;
			}
			else {
				H1(i, j) = dLnphiLdLnK - dLnphiGdLnK;
			}
		}
	}

	for (i = 0; i < Nc; i++) {
		H1(Nc, i) = ZI[i] / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		//H1(Nc, i) = 0.0;
	}
	df = 0.0;
	for (i = 0; i < Nc; i++) {
		df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
	}
	H1(Nc, Nc) = df;
	H0 = H1.inverse();
	MatrixXd ttt = H0 * H1;
	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		Lnk0 = LnK1;

		LnK1 = Lnk0 - 0.5 * H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		nv = LnK1(Nc);
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);
		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}
		for (i = 0; i < Nc; i++) {
			g1(i) = LnK1(i) + LnphiL(i) - LnphiG(i);
		}
		f = 0.0;
		for (i = 0; i < Nc; i++) {
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
		}
		LnK1(Nc) = nv;
		g1(Nc) = f;


		s1 = LnK1 - Lnk0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc + 1);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc + 1; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc + 1; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc + 1; i++) {
			for (j = 0; j < Nc + 1; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 0.999999 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.000001 / (1.0 - Kmax);
		}



		/*	if (nv > 1.0) {
				nv = 0.999999;
			}
			if (nv < 0) {
				nv = 0.000001;
			}*/

		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
	}

	return 0;
}





void Compositional::zeroInitializatinoCubicEOS(double a, double b, double* S, double A, double B, double Z, string phase) {
	register int i;
	phase = "";
	a = 0.0;
	b = 0.0;
	A = 0.0;
	B = 0.0;
	for (i = 0; i < Nc; i++) {
		S[i] = 0.0;
	}
	Z = 0.0;

}


double Compositional::twoPhaseStabilityFlash(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiZ = new double[Nc];
	double* phiU = new double[Nc];
	double* fZ = new double[Nc];
	double* fU = new double[Nc];
	double* u = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiZ[i] = 0.0;
		phiU[i] = 0.0;
		S[i] = 0.0;
		fZ[i] = 0.0;
		fU[i] = 0.0;
		u[i] = 0.0;
	}
	VectorXd g1(Nc), g0(Nc), LnK1(Nc), Lnk0(Nc), LnphiZ(Nc),
		s1(Nc), y1(Nc), LnphiU(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);

	for (i = 0; i < Nc; i++) {
		u[i] = K[i] * ZI[i];
	}
	normalizArr(u, Nc);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiZ(i) = log(phiZ[i]);
		LnphiU(i) = log(phiU[i]);
		LnK1(i) = log(K[i]);
	}

	g1 = LnK1 + LnphiU - LnphiZ;

	g0 = g1;
	Lnk0 = LnK1;
	LnK1 = Lnk0 - H0 * g0;
	double sumEst = LnK1.array().abs().sum();
	if (AlmostEqual(sumEst, 0, Nc * 1e-1)) {
		for (i = 0; i < Nc; i++) {
			u[i] = 1.0 / K[i] * ZI[i];
		}
		normalizArr(u, Nc);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

		H0.setIdentity(Nc, Nc);
		for (i = 0; i < Nc; i++) {
			LnphiZ(i) = log(phiZ[i]);
			LnphiU(i) = log(phiU[i]);
			LnK1(i) = log(K[i]);
		}

		g1 = LnK1 + LnphiU - LnphiZ;

		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - H0 * g0;
	}


	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;

		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		for (i = 0; i < Nc; i++) {
			u[i] = K[i] * ZI[i];
		}
		normalizArr(u, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);
		for (i = 0; i < Nc; i++) {
			LnphiU(i) = log(phiU[i]);
			LnphiZ(i) = log(phiZ[i]);
		}

		g1 = LnK1 + LnphiU - LnphiZ;

		s1 = LnK1 - Lnk0;

		y1 = g1 - g0;
		H1 = H0 + (s1 - H0 * y1) * (s1.transpose() * H0);
		MatrixXd temp = (s1.transpose() * H0 * y1);

		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				H1(i, j) = H1(i, j) / temp(0, 0);
			}
		}

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= f / df;
		err_2 = abs(nv - nv_old);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 0.98 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.02 / (1.0 - Kmax);
		}

		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - H0 * g0;
	}
	double Su = 0;
	for (i = 0; i < Nc; i++) {
		K[i] = exp(LnK1(i));
		Su += K[i] * ZI[i];
	}
	return Su;


}


double Compositional::twoPhaseStabilityFlash_BFGS(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiZ = new double[Nc];
	double* phiU = new double[Nc];
	double* fZ = new double[Nc];
	double* fU = new double[Nc];
	double* u = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 0; i < Nc; i++) {
		phiZ[i] = 0.0;
		phiU[i] = 0.0;
		S[i] = 0.0;
		fZ[i] = 0.0;
		fU[i] = 0.0;
		u[i] = 0.0;
		K_old[i] = K[i];
	}
	VectorXd g1(Nc), g0(Nc), LnK1(Nc), Lnk0(Nc), LnphiZ(Nc),
		s1(Nc), y1(Nc), LnphiU(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);

	for (i = 0; i < Nc; i++) {
		u[i] = K[i] * ZI[i];
	}
	normalizArr(u, Nc);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiZ(i) = log(phiZ[i]);
		LnphiU(i) = log(phiU[i]);
		LnK1(i) = log(K[i]);
	}

	g1 = LnK1 + LnphiU - LnphiZ;

	g0 = g1;
	Lnk0 = LnK1;
	LnK1 = Lnk0 - H0 * g0;
	double sumEst = LnK1.array().abs().sum();
	if (AlmostEqual(sumEst, 0, Nc * 2 * 1e-1)) {
		for (i = 0; i < Nc; i++) {
			K[i] = 1 / K_old[i];
			u[i] = K[i] * ZI[i];
		}
		normalizArr(u, Nc);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

		H0.setIdentity(Nc, Nc);
		for (i = 0; i < Nc; i++) {
			LnphiZ(i) = log(phiZ[i]);
			LnphiU(i) = log(phiU[i]);
			LnK1(i) = log(K[i]);
		}

		g1 = LnK1 + LnphiU - LnphiZ;

		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - H0 * g0;
	}
	//Kmin = findMin(K, Nc);
	//Kmax = findMax(K, Nc);
	//fstream testnv;
	//testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
	//nv = -1.0;
	//int itrtest = 0;
	//while (nv < Kmax + 1.0) {
	//	itrtest++;
	//	f = 0.0;
	//	df = 0.0;
	//	for (i = 0; i < Nc; i++) {
	//		f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
	//		df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
	//	}
	//	testnv << itrtest << "		" << nv << "		" << f << "			" << df << "			" << endl;
	//	nv = nv + 0.002;
	//}


	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;

		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		for (i = 0; i < Nc; i++) {
			u[i] = K[i] * ZI[i];
		}
		normalizArr(u, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);
		for (i = 0; i < Nc; i++) {
			LnphiU(i) = log(phiU[i]);
			LnphiZ(i) = log(phiZ[i]);
		}

		g1 = LnK1 + LnphiU - LnphiZ;

		s1 = LnK1 - Lnk0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;



		nv_old = nv;
		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);
		/*	fstream testnv;
			testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
			nv = Kmin;
			int itrtest = 0;
			while (nv < Kmax) {
				itrtest++;
				f = 0.0;
				df = 0.0;
				for (i = 0; i < Nc; i++) {
					f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
					df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
				}
				testnv << itrtest << "		" << nv << "		" << f << "			" << df << "			" << endl;
				nv = nv + 0.002;
			}*/


		err_2 = 1.0;
		int itr_2 = 0;
		/*while ((err_2 > 1e-8) && (itr_2 < 100)){
			itr_2++;*/


			/*fstream testnv;
			testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
			nv = 0.0;
			int itrtest = 0;
			while (nv < Kmax) {
				itrtest++;
				f = 0.0;
				df = 0.0;
				for (i = 0; i < Nc; i++) {
					f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
					df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
				}
				testnv << itrtest << "		" << nv << "		" << f << "			" << df << "			" << endl;
				nv = nv + 0.002;
			}*/



		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
			df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
		}
		nv_old = nv;
		nv -= 0.01 * f / df;
		//err_2 = abs(nv - nv_old);




		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 0.98 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.02 / (1.0 - Kmax);
		}
		//}
		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - H0 * g0;
	}
	double Su = 0;
	for (i = 0; i < Nc; i++) {
		K[i] = exp(LnK1(i));
		Su += K[i] * ZI[i];
	}
	return Su;


}


double Compositional::twoPhaseStabilityFlash_BFGS_FullyImplicit(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double* K, double& nv) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiZ = new double[Nc];
	double* phiU = new double[Nc];
	double* fZ = new double[Nc];
	double* fU = new double[Nc];
	double* u = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 0; i < Nc; i++) {
		phiZ[i] = 0.0;
		phiU[i] = 0.0;
		S[i] = 0.0;
		fZ[i] = 0.0;
		fU[i] = 0.0;
		u[i] = 0.0;
		K_old[i] = K[i];
	}
	VectorXd g1(Nc + 1), g0(Nc + 1), LnK1(Nc + 1), Lnk0(Nc + 1), LnphiZ(Nc),
		s1(Nc + 1), y1(Nc + 1), LnphiU(Nc);
	MatrixXd H1(Nc + 1, Nc + 1), H0(Nc + 1, Nc + 1);

	for (i = 0; i < Nc; i++) {
		u[i] = K[i] * ZI[i];
	}
	normalizArr(u, Nc);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

	H0.setIdentity(Nc + 1, Nc + 1);
	for (i = 0; i < Nc; i++) {
		LnphiZ(i) = log(phiZ[i]);
		LnphiU(i) = log(phiU[i]);
		LnK1(i) = log(K[i]);
	}
	for (i = 0; i < Nc; i++) {
		g1(i) = LnK1(i) + LnphiU(i) - LnphiZ(i);
	}

	f = 0.0;
	for (i = 0; i < Nc; i++) {
		f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
	}
	LnK1(Nc) = nv;
	g1(Nc) = f;



	g0 = g1;
	Lnk0 = LnK1;
	LnK1 = Lnk0 - H0 * g0;
	double sumEst = LnK1.array().abs().sum();
	if (AlmostEqual(sumEst, 0, Nc * 2 * 1e-1)) {
		for (i = 0; i < Nc; i++) {
			K[i] = 1 / K_old[i];
			u[i] = K[i] * ZI[i];
		}
		normalizArr(u, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);

		H0.setIdentity(Nc, Nc);
		for (i = 0; i < Nc; i++) {
			LnphiZ(i) = log(phiZ[i]);
			LnphiU(i) = log(phiU[i]);
			LnK1(i) = log(K[i]);
		}
		f = 0.0;
		for (i = 0; i < Nc; i++) {
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
		}
		LnK1(Nc) = nv;
		g1(Nc) = f;
		for (i = 0; i < Nc; i++) {
			g1(i) = LnK1(i) + LnphiU(i) - LnphiZ(i);
		}


		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - H0 * g0;
	}



	//Kmin = findMin(K, Nc);
	//Kmax = findMax(K, Nc);
	//fstream testnv;
	//testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
	//nv = -1.0;
	//int itrtest = 0;
	//while (nv < Kmax + 1.0) {
	//	itrtest++;
	//	f = 0.0;
	//	df = 0.0;
	//	for (i = 0; i < Nc; i++) {
	//		f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
	//		df += -(ZI[i] * (K[i] - 1.0) * (K[i] - 1.0)) / ((1.0 + (K[i] - 1.0) * nv) * (1.0 + (K[i] - 1.0) * nv));
	//	}
	//	testnv << itrtest << "		" << nv << "		" << f << "			" << df << "			" << endl;
	//	nv = nv + 0.002;
	//}


	err_1 = 1.0;
	int itr = 0;
	while ((err_1 > 1e-18) && (itr < 10000)) {
		itr++;

		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnK1(i));
		}
		for (i = 0; i < Nc; i++) {
			u[i] = K[i] * ZI[i];
		}
		nv = LnK1(Nc);
		normalizArr(u, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, ZI, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiZ);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, u, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiU);
		for (i = 0; i < Nc; i++) {
			LnphiU(i) = log(phiU[i]);
			LnphiZ(i) = log(phiZ[i]);
		}
		f = 0.0;
		for (i = 0; i < Nc; i++) {
			f += (ZI[i] * (K[i] - 1.0)) / (1.0 + (K[i] - 1.0) * nv);
		}
		LnK1(Nc) = nv;
		g1(Nc) = f;

		for (i = 0; i < Nc; i++) {
			g1(i) = LnK1(i) + LnphiU(i) - LnphiZ(i);
		}


		s1 = LnK1 - Lnk0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc + 1);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc + 1; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc + 1; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc + 1; i++) {
			for (j = 0; j < Nc + 1; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);

		if (nv > 1.0 / (1.0 - Kmin)) {
			nv = 0.999999 / (1.0 - Kmin);
		}
		if (nv < 1.0 / (1.0 - Kmax)) {
			nv = 1.000001 / (1.0 - Kmax);
		}
		//}
		auto tempvec = (Lnk0 - LnK1).array() / Lnk0.array();
		err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		g0 = g1;
		Lnk0 = LnK1;
		LnK1 = Lnk0 - 0.5 * H0 * g0;
	}
	double Su = 0;
	for (i = 0; i < Nc; i++) {
		K[i] = exp(LnK1(i));
		Su += K[i] * ZI[i];
	}
	return Su;


}




double Compositional::enthalpyReal(double* x, double T, double P, double* a0, double* a1, double* a2, double* a3, double* a4,
	double Z, double A, double B, double phi) {
	double HReal;
	double HIdeal;
	double HRes;
	HIdeal = enthalpyIdeal(x, T, a0, a1, a2, a3, a4);
	double dZdT = dZdTFunc(Z, A, B, T);
	double dphidT = dphidTFunc(dZdT, A, B, Z, T);
	HRes = enthalpyResedual(T, phi, dphidT);
	HReal = HIdeal + HRes;
	return HReal;
}

double Compositional::enthalpyIdeal(double* x, double T, double* a0, double* a1, double* a2, double* a3, double* a4) {
	double HIdeal = 0.0;
	register int i;

	for (i = 0; i < Nc; i++) {
		HIdeal += ((a0[i] * T + a1[i] * T * T / 2.0 * 1e-3 + a2[i] * T * T * T / 3.0 * 1e-5 +
			a3[i] * T * T * T * T / 4.0 * 1e-8 + a4[i] * T * T * T * T * T / 5.0 * 1e-11) -
			(a0[i] * Tref_ + a1[i] * Tref_ * Tref_ / 2.0 * 1e-3 + a2[i] * Tref_ * Tref_ * Tref_ / 3.0 * 1e-5 +
				a3[i] * Tref_ * Tref_ * Tref_ * Tref_ / 4.0 * 1e-8 + a4[i] * Tref_ * Tref_ * Tref_ * Tref_ * Tref_ / 5.0 * 1e-11))* x[i];
	}

	return Rgas * HIdeal;
}


double Compositional::HeatCapacityIdeal(double* x, double T, double* a0, double* a1, double* a2, double* a3, double* a4) {
	double cpIdeal = 0.0;
	register int i;
	for (i = 0; i < Nc; i++) {
		cpIdeal += (a0[i] + a1[i] * T * 1e-3 + a2[i] * T * T * 1e-5 + a3[i] * T * T * T * 1e-8 + a4[i] * T * T * T * T * 1e-11) * x[i];
	}


	return Rgas * cpIdeal;
}

double Compositional::dZdTFunc(double Z, double A, double B, double T) {
	double da1dT, da2dT, da3dT;
	double c;
	if ((EOS == "RK") || (EOS == "SRK") || (EOS == "GD")) {
		c = 0.0;
	}
	else if ((EOS == "PR") || ("PR78")) {
		c = 1.0;
	}
	double dAdT = A * (-2.0 / T);
	double dBdT = B / (-T);
	double a1 = -(1.0 - c * B);
	double a2 = A - B * (1.0 + c) - B * B * (1.0 + 2.0 * c);
	da1dT = c * dBdT;
	da2dT = dAdT - dBdT * (1.0 + c) - 2.0 * B * dBdT * (1.0 + 2.0 * c);
	da3dT = -(B * dAdT + A * dBdT - c * (3.0 * B * B * dBdT + 2.0 * B * dBdT));
	return (-Z * Z * da1dT - Z * da2dT - da3dT) / (3 * Z * Z + 2 * Z * a1 + a2);
}

double Compositional::dphidTFunc(double dZdT, double A, double B, double Z, double T) {
	double dphi;
	double d1, d2;
	if ((EOS == "RK") || (EOS == "SRK") || (EOS == "GD")) {
		d1 = 0.0;
		d2 = 1.0;
	}
	else if ((EOS == "PR") || ("PR78")) {
		d1 = 1 - sqrt(2);
		d2 = 1 + sqrt(2);
	}
	double dAdT = A * (-2.0 / T);
	double dBdT = B / (-T);
	dphi = dZdT - dZdT / (Z - B) - 2.0 / (d2 - d1) * (dAdT / B - A / (B * B) * dBdT) * log((Z + d2 * B) / (Z + d1 * B)) -
		2.0 / (d2 - d1) * (A / B) *
		(((dZdT + d2 * dBdT) * (Z + d1 * B) - (dZdT + d1 * dBdT) * (Z + d2 * B)) / ((Z + d2 * B) * (Z + d1 * B)));
	return dphi;
}

double Compositional::enthalpyResedual(double T, double phi, double dphidT) {

	double HRes = -Rgas * T * T * dphidT / phi;
	return HRes;
}

void Compositional::setNc(int N) {
	Nc = N;
}
void Compositional::setEOS(string inEOS) {
	EOS = inEOS;
}


double Compositional::viscosity_LBC(double P, double T, double rho, double* x, double* Tc, double* Pc, double* Vc, double* Mw) {
	double* etaStar_i = new double[Nc];
	double* zeta_i = new double[Nc];
	double* Tr = new double[Nc];
	double a1 = 0.10230;
	double a2 = 0.023364;
	double a3 = 0.058533;
	double a4 = -0.040758;
	double a5 = 0.0093324;
	register int i;

	double rhoC = 0.0;
	double TcMix = 0.0;
	double MwMix = 0.0;
	double PcMix = 0.0;
	for (i = 0; i < Nc; i++) {
		rhoC += x[i] * Vc[i];
		TcMix += x[i] * Tc[i];
		MwMix += x[i] * Mw[i];
		PcMix += x[i] * Pc[i];
		zeta_i[i] = pow(Tc[i], 1.0 / 6.0) / sqrt(Mw[i] * 1000.0) / pow(cbrtl(Pas_to_atm(Pc[i])), 2.0);
		Tr[i] = T / Tc[i];
	}
	rhoC = 1 / rhoC;
	double zeta = pow(TcMix, 1.0 / 6.0) / sqrt(MwMix * 1000.0) / pow(cbrtl(Pas_to_atm(PcMix)), 2.0);



	for (i = 0; i < Nc; i++) {
		if (Tr[i] <= 1.5) {
			etaStar_i[i] = 34 * 1e-5 / zeta_i[i] * pow(Tr[i], 0.94);
		}
		else if (Tr[i] > 1.5) {
			etaStar_i[i] = 17.78 * 1e-5 / zeta_i[i] * pow(4.58 * Tr[i] - 1.67, 5.0 / 8.0);
		}
	}


	double sum1 = 0.0;
	double sum2 = 0.0;
	for (i = 0; i < Nc; i++) {
		sum1 += x[i] * etaStar_i[i] * sqrt(Mw[i]);
		sum2 += x[i] * sqrt(Mw[i]);
	}
	double etaStar = sum1 / sum2;


	double rhoR = rho / MwMix / rhoC;
	double eta = (pow((a1 + a2 * rhoR + a3 * rhoR * rhoR + a4 * rhoR * rhoR * rhoR + a5 * rhoR * rhoR * rhoR * rhoR), 4.0) - 1e-4) / zeta + etaStar;

	return eta;
}

double Compositional::density_f(double P, double T, double MwMix, double Z) {

	double d = P * MwMix / (Z * Rgas * T);
	return d;
}

double Compositional::molecularWight(double* x, double* Mw) {
	register int i;
	double MwMix = 0.0;

	for (i = 0; i < Nc; i++)
	{
		MwMix += x[i] * Mw[i];
	}
	return MwMix;
}


double Compositional::GasOilIFT(double rhoL, double rhoG, double* Par, double* x, double* y) {
	double IFT = 0.0;
	double sum = 0.0;
	register int i;
	for (i = 0; i < Nc; i++) {
		sum += rhoL * Par[i] * x[i] - rhoG * Par[i] * y[i];
	}
	IFT = sum * sum * sum * sum;

	return IFT;
}


int Compositional::twoPhaseFlashEnv_BFGS_MSS_bubble(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* phiL_f = new double[Nc];
	double* phiG_f = new double[Nc];
	double* phiL_b = new double[Nc];
	double* phiG_b = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double* K = new double[Nc];
	double* dK = new double[Nc];
	double* dK_f = new double[Nc];
	double* dK_b = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc), g0(Nc), LnAlpha1(Nc), LnAlpha0(Nc), LnphiL(Nc),
		s1(Nc), s1t(Nc), y1(Nc), LnphiG(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);
	wilsonK(P, T, Tc, Pc, ACF, K);
	twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, ZI, K, nv);
	nv = 0.0;
	double errT = 1.0;
	//double T_old;
	double P_old;
	//while (errT > 1e-8) {
	//	wilsonK(P, T, Tc, Pc, ACF, K);
	//	//dTwilsonK(P, T, Tc, Pc, ACF, dK);
	//	dPwilsonK(P, T, Tc, Pc, ACF, dK);
	//	f = 0.0;
	//	df = 0.0;
	//	inner_counter = 0;
	//	for (i = 0; i < Nc; i++) {
	//		inner_counter++;
	//		f += ZI[i] * (K[i] - 1.0);
	//		df += ZI[i] * dK[i];
	//	}
	//	/*T_old = T;
	//	T -= f / df;
	//	errT = abs((T - T_old) / T_old);*/
	//	P_old = P;
	//	P -= f / df;
	//	errT = abs((P - P_old) / P_old);
	//}

	//g1(Nc) = 0.0;
	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
		//g1(Nc) += y[i] - x[i];
	}
	//g1(Nc) = log(g1(Nc));

	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnAlpha1(i) = log(K[i]);
		g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
	}
	//LnAlpha1(Nc) = log(T);
	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		LnAlpha0 = LnAlpha1;

		LnAlpha1 = LnAlpha0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnAlpha1(i));
		}
		//T = exp(LnAlpha1(Nc));
		//g1(Nc) = 0.0;
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
			//g1(Nc) += y[i] - x[i];
		}

		//g1(Nc) = log(g1(Nc));
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
			g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
		}

		s1 = LnAlpha1 - LnAlpha0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;



		/*for (i = 0; i < Nc; i++) {
			K[i] = exp(LnphiL[i] - LnphiG[i]);
		}*/

		P_old = P;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*1.001, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*0.999, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*1.001, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*0.999, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);


		/*T_old = T;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);
*/
/*for (i = 0; i < Nc; i++) {
	dK[i] = (phiL_f[i] / phiG_f[i] - phiL_b[i] / phiG_b[i]) / (2 * 0.001 * T);
}*/

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);
		//fstream testnv;
		//testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
		//P = 1e+6;
		//int itrtest = 0;
		//while (P < 3 * 1e+6) {
		//	itrtest++;
		//	f = 0.0;
		//	df = 0.0;
		//	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		//	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		//	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		//	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		//	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		//	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		//	for (i = 0; i < Nc; i++) {
		//		LnphiG(i) = log(phiG[i]);
		//		LnphiL(i) = log(phiL[i]);
		//		K[i] = phiL[i] / phiG[i];
		//	}
		//	for (i = 0; i < Nc; i++) {
		//		f += ZI[i] * (K[i] - 1.0);
		//		//df += ZI[i] * dK[i];
		//		df += ZI[i] * K[i] * ((phiL_f[i] - phiL_b[i]) / (2 * 0.001 * P) - (phiG_f[i] - phiG_b[i]) / (2 * 0.001 * P));
		//	}
		//	testnv << itrtest << "		" << P << "		" << f << "			" << df << "			" << endl;
		//	P = P + 1e+4;
		//}
		//
		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += ZI[i] * (K[i] - 1.0);
			//df += ZI[i] * dK[i];
			df += ZI[i] * K[i] * ((phiL_f[i] - phiL_b[i]) / (2 * 0.001 * P) - (phiG_f[i] - phiG_b[i]) / (2 * 0.001 * P));
		}
		P_old = P;
		P -= f / df;

		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}



		//auto tempvec = (LnAlpha0 - LnAlpha1).array() / LnAlpha0.array();
		//err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		//err_1 = abs((T_old - T) / T_old);

		err_1 = abs((P_old - P) / P_old);
	}

	return 0;
}

int Compositional::twoPhaseFlashEnv_BFGS_MSS_bubble2(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* phiL_f = new double[Nc];
	double* phiG_f = new double[Nc];
	double* phiL_b = new double[Nc];
	double* phiG_b = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double* K = new double[Nc];
	double* dK = new double[Nc];
	double* dK_f = new double[Nc];
	double* dK_b = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc + 1), g0(Nc + 1), LnAlpha1(Nc + 1), LnAlpha0(Nc + 1), LnphiL(Nc),
		s1(Nc + 1), s1t(Nc + 1), y1(Nc + 1), LnphiG(Nc);
	MatrixXd H1(Nc + 1, Nc + 1), H0(Nc + 1, Nc + 1);
	wilsonK(P, T, Tc, Pc, ACF, K);
	twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, ZI, K, nv);
	nv = 0.0;
	double errT = 1.0;
	//double T_old;
	double P_old;
	//while (errT > 1e-8) {
	//	wilsonK(P, T, Tc, Pc, ACF, K);
	//	//dTwilsonK(P, T, Tc, Pc, ACF, dK);
	//	dPwilsonK(P, T, Tc, Pc, ACF, dK);
	//	f = 0.0;
	//	df = 0.0;
	//	inner_counter = 0;
	//	for (i = 0; i < Nc; i++) {
	//		inner_counter++;
	//		f += ZI[i] * (K[i] - 1.0);
	//		df += ZI[i] * dK[i];
	//	}
	//	/*T_old = T;
	//	T -= f / df;
	//	errT = abs((T - T_old) / T_old);*/
	//	P_old = P;
	//	P -= f / df;
	//	errT = abs((P - P_old) / P_old);
	//}

	g1(Nc) = 0.0;
	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
		g1(Nc) += y[i] - x[i];
	}
	//g1(Nc) = log(g1(Nc));

	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc + 1, Nc + 1);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnAlpha1(i) = log(K[i]);
		g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
	}
	LnAlpha1(Nc) = log(P);
	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		LnAlpha0 = LnAlpha1;

		LnAlpha1 = LnAlpha0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnAlpha1(i));
		}
		P_old = P;
		P = exp(LnAlpha1(Nc));
		g1(Nc) = 0.0;
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
			g1(Nc) += y[i] - x[i];
		}

		//g1(Nc) = log(g1(Nc));
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
			g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
		}

		s1 = LnAlpha1 - LnAlpha0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc + 1);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc + 1; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc + 1; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc + 1; i++) {
			for (j = 0; j < Nc + 1; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;



		/*for (i = 0; i < Nc; i++) {
			K[i] = exp(LnphiL[i] - LnphiG[i]);
		}*/


		/*	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P*1.001, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P*0.999, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P*1.001, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T, P*0.999, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);*/


			/*T_old = T;


			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T*1.001, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T*0.999, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T*1.001, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

			zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
			cubicEOS(T*0.999, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
			cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);
	*/
	/*for (i = 0; i < Nc; i++) {
		dK[i] = (phiL_f[i] / phiG_f[i] - phiL_b[i] / phiG_b[i]) / (2 * 0.001 * T);
	}*/

	////f = 0.0;
	////df = 0.0;
	////inner_counter = 0;
	////for (i = 0; i < Nc; i++) {
	////	inner_counter++;
	////	f += ZI[i] * (K[i] - 1.0);
	////	//df += ZI[i] * dK[i];
	////	df += ZI[i] * K[i] * ((phiL_f[i] - phiL_b[i]) / (2 * 0.001 * P) - (phiG_f[i] - phiG_b[i]) / (2 * 0.001 * P));
	////}
	////P_old = P;
	////P -= f / df;

	/*for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
	}
	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
	}*/



	//auto tempvec = (LnAlpha0 - LnAlpha1).array() / LnAlpha0.array();
	//err_1 = tempvec.array().square().sum();
	//err_1 = tempvec.no
	//err_1 = abs((T_old - T) / T_old);

		err_1 = abs((P_old - P) / P_old);
	}

	return 0;
}




int Compositional::twoPhaseFlashEnv_BFGS_MSS_dew(double P, double T, double* ai, double* bi, double** BIC, double* ZI,
	double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* phiL_f = new double[Nc];
	double* phiG_f = new double[Nc];
	double* phiL_b = new double[Nc];
	double* phiG_b = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double* K = new double[Nc];
	double* dK = new double[Nc];
	double* dK_f = new double[Nc];
	double* dK_b = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc), g0(Nc), LnAlpha1(Nc), LnAlpha0(Nc), LnphiL(Nc),
		s1(Nc), s1t(Nc), y1(Nc), LnphiG(Nc);
	MatrixXd H1(Nc, Nc), H0(Nc, Nc);
	wilsonK(P, T, Tc, Pc, ACF, K);
	twoPhaseStabilityFlash_BFGS_FullyImplicit(P, T, ai, bi, BIC, ZI, K, nv);
	nv = 1.0;
	double errT = 1.0;
	double T_old;
	//double P_old;
	//while (errT > 1e-8) {
	//	wilsonK(P, T, Tc, Pc, ACF, K);
	//	dTwilsonK(P, T, Tc, Pc, ACF, dK);
	//	//dPwilsonK(P, T, Tc, Pc, ACF, dK);
	//	f = 0.0;
	//	df = 0.0;
	//	inner_counter = 0;
	//	for (i = 0; i < Nc; i++) {
	//		inner_counter++;
	//		f += ZI[i] / K[i] - y[i];
	//		df += - ZI[i] * dK[i] / (K[i] * K[i]);
	//	}
	//	T_old = T;
	//	T -= f / df;
	//	errT = abs((T - T_old) / T_old);
	//	/*P_old = P;
	//	P -= f / df;
	//	errT = abs((P - P_old) / P_old);*/
	//}

	//g1(Nc) = 0.0;
	for (i = 0; i < Nc; i++) {
		y[i] = ZI[i];
		x[i] = y[i] / K[i];
		//g1(Nc) += y[i] - x[i];
	}
	//g1(Nc) = log(g1(Nc));

	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc, Nc);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnAlpha1(i) = log(K[i]);
		g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
	}
	//LnAlpha1(Nc) = log(T);
	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		LnAlpha0 = LnAlpha1;

		LnAlpha1 = LnAlpha0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnAlpha1(i));
		}
		//T = exp(LnAlpha1(Nc));
		//g1(Nc) = 0.0;
		for (i = 0; i < Nc; i++) {
			y[i] = ZI[i];
			x[i] = y[i] / K[i];
			//g1(Nc) += y[i] - x[i];
		}

		//g1(Nc) = log(g1(Nc));
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
			g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
		}

		s1 = LnAlpha1 - LnAlpha0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc; i++) {
			for (j = 0; j < Nc; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;



		/*for (i = 0; i < Nc; i++) {
			K[i] = exp(LnphiL[i] - LnphiG[i]);
		}*/

		T_old = T;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);


		/*T_old = T;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);
*/
/*for (i = 0; i < Nc; i++) {
	dK[i] = (phiL_f[i] / phiG_f[i] - phiL_b[i] / phiG_b[i]) / (2 * 0.001 * T);
}*/

		Kmin = findMin(K, Nc);
		Kmax = findMax(K, Nc);
		//fstream testnv;
		//testnv.open("D:/well_res/coupled-well-res/WellboreReservoir_CPP/testnv.txt", ios::out);
		//T = 200.0;
		//int itrtest = 0;
		//while (T < 700.0) {
		//	itrtest++;
		//	f = 0.0;
		//	df = 0.0;
		//	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		//	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		//	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		//	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		//	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		//	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		//	for (i = 0; i < Nc; i++) {
		//		LnphiG(i) = log(phiG[i]);
		//		LnphiL(i) = log(phiL[i]);
		//		K[i] = phiL[i] / phiG[i];
		//	}
		//	for (i = 0; i < Nc; i++) {
		//		f += ZI[i] / K[i] - 1.0;
		//		//df += ZI[i] * dK[i];
		//		df += ZI[i] / K[i] * ((log(phiG_f[i]) - log(phiG_b[i])) / (2 * 0.001 * T) - (log(phiL_f[i]) - log(phiL_b[i])) / (2 * 0.001 * T));
		//	}
		//	testnv << itrtest << "		" << T << "		" << f << "			" << df << "			" << endl;
		//	T = T + 1.0;
		//}



		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += ZI[i] / K[i] - 1.0;
			//df += ZI[i] * dK[i];
			df += ZI[i] / K[i] * ((log(phiG_f[i]) - log(phiG_b[i])) / (2 * 0.001 * T) - (log(phiL_f[i]) - log(phiL_b[i])) / (2 * 0.001 * T));
		}
		/*P_old = P;
		P -= f / df;*/
		T_old = T;
		T -= f / df;
		//errT = abs((T - T_old) / T_old);

		for (i = 0; i < Nc; i++) {
			y[i] = ZI[i];
			x[i] = y[i] / K[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}



		//auto tempvec = (LnAlpha0 - LnAlpha1).array() / LnAlpha0.array();
		//err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		//err_1 = abs((T_old - T) / T_old);
		err_1 = abs((T_old - T) / T_old);
	}

	return 0;
}






int Compositional::twoPhaseFlashEnv_BFGS_MSS2(double P, double P_old0, double T, double T_old, double* ai, double* bi, double** BIC, double* ZI,
	double& nv, double* x, double* y, double* Tc, double* Pc, double* ACF) {
	register int i, j;
	double tol_1 = 1e-18;
	double tol_2 = 1e-8;
	double err_1 = 1.0;
	double err_2 = 1.0;
	double nv_old;
	double* K_old = new double[Nc];
	double* S = new double[Nc];
	double* phiL = new double[Nc];
	double* phiG = new double[Nc];
	double* phiL_f = new double[Nc];
	double* phiG_f = new double[Nc];
	double* phiL_b = new double[Nc];
	double* phiG_b = new double[Nc];
	double* fL = new double[Nc];
	double* fG = new double[Nc];
	double* K = new double[Nc];
	double* dK = new double[Nc];
	double* dK_f = new double[Nc];
	double* dK_b = new double[Nc];
	double f, df, inner_counter;
	double a = 0;
	double b = 0.0;
	double A = 0.0;
	double B = 0.0;
	double Z = 0.0;
	string phase;
	double Kmax, Kmin;
	for (i = 1; i < Nc; i++) {
		phiL[i] = 0.0;
		phiG[i] = 0.0;
		S[i] = 0.0;
		fL[i] = 0.0;
		fG[i] = 0.0;
	}
	VectorXd g1(Nc + 1), g0(Nc + 1), LnAlpha1(Nc + 1), LnAlpha0(Nc + 1), LnphiL(Nc),
		s1(Nc + 1), s1t(Nc + 1), y1(Nc + 1), LnphiG(Nc);
	MatrixXd H1(Nc + 1, Nc + 1), H0(Nc + 1, Nc + 1);
	wilsonK(P, T, Tc, Pc, ACF, K);
	twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, ZI, K, nv);
	nv = 0;
	double errT = 1.0;
	//double T_old;
	double P_old = P;
	//double P_old0;
	//while (errT > 1e-8) {
	//	wilsonK(P, T, Tc, Pc, ACF, K);
	//	//dTwilsonK(P, T, Tc, Pc, ACF, dK);
	//	dPwilsonK(P, T, Tc, Pc, ACF, dK);
	//	f = 0.0;
	//	df = 0.0;
	//	inner_counter = 0;
	//	for (i = 0; i < Nc; i++) {
	//		inner_counter++;
	//		f += ZI[i] * (K[i] - 1.0);
	//		df += ZI[i] * dK[i];
	//	}
	//	/*T_old = T;
	//	T -= f / df;
	//	errT = abs((T - T_old) / T_old);*/
	//	P_old = P;
	//	P -= f / df;
	//	errT = abs((P - P_old) / P_old);
	//}

	//g1(Nc) = 0.0;
	for (i = 0; i < Nc; i++) {
		x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
		y[i] = K[i] * x[i];
		//g1(Nc) += y[i] - x[i];
	}
	//g1(Nc) = log(g1(Nc));

	normalizArr(x, Nc);
	normalizArr(y, Nc);
	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

	zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
	cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
	cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

	H0.setIdentity(Nc + 1, Nc + 1);
	for (i = 0; i < Nc; i++) {
		LnphiG(i) = log(phiG[i]);
		LnphiL(i) = log(phiL[i]);
		LnAlpha1(i) = log(K[i]);
		g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
	}
	LnAlpha1(Nc) = log(P);
	g1(Nc) = ((P_old0 - P) / 1e+6) * ((P_old0 - P) / 1e+6) + ((T_old - T) / 1e+4) * ((T_old - T) / 1e+4) - 1;
	err_1 = 1.0;
	int itr = 0;
	while (err_1 > tol_1) {
		itr++;
		g0 = g1;
		LnAlpha0 = LnAlpha1;

		LnAlpha1 = LnAlpha0 - H0 * g0;
		for (i = 0; i < Nc; i++) {
			K[i] = exp(LnAlpha1(i));
		}
		T = exp(LnAlpha1(Nc));
		//g1(Nc) = 0.0;
		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
			//g1(Nc) += y[i] - x[i];
		}

		//g1(Nc) = log(g1(Nc));
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
			g1(i) = LnAlpha1(i) - LnphiL(i) + LnphiG(i);
		}
		g1(Nc) = log10(abs((P_old0 - P) / 1e+5)) * log10(abs((P_old0 - P) / 1e+5)) + log10(abs(T_old - T)) * log10(abs(T_old - T)) - 1;
		s1 = LnAlpha1 - LnAlpha0;
		y1 = g1 - g0;
		MatrixXd temp1 = s1.transpose()*y1; // s1'*y1
		MatrixXd temp2 = y1.transpose() * H0; // *y1;		// y1'*H0*y1
		VectorXd temp3 = H0 * y1;
		VectorXd BFGS_U(Nc + 1);
		MatrixXd s1Ps1t = s1 * s1.transpose();	//s1*s1'
		MatrixXd y1Py1t = y1 * y1.transpose();	//y1*y1'
		MatrixXd H0Py1Py1tPH0 = H0 * y1Py1t * H0; // H0*y1*y1'*H0
		for (i = 0; i < Nc + 1; i++) {
			BFGS_U(i) = s1(i) / temp1(0, 0) - temp3(i) / temp2(0, 0);
			for (j = 0; j < Nc + 1; j++) {
				s1Ps1t(i, j) = s1Ps1t(i, j) / temp1(0, 0);
				H0Py1Py1tPH0(i, j) = H0Py1Py1tPH0(i, j) / temp2(0, 0);
			}
		}

		MatrixXd BFGS_UU = BFGS_U * BFGS_U.transpose();
		for (i = 0; i < Nc + 1; i++) {
			for (j = 0; j < Nc + 1; j++) {
				BFGS_UU(i, j) = BFGS_UU(i, j) * temp2(0, 0);
			}
		}

		H1 = H0 + s1Ps1t - H0Py1Py1tPH0 + BFGS_UU;



		/*for (i = 0; i < Nc; i++) {
			K[i] = exp(LnphiL[i] - LnphiG[i]);
		}*/

		P_old0 = P;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*1.001, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*0.999, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*1.001, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P*0.999, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);


		/*T_old = T;


		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL_b);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*1.001, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_f);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T*0.999, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG_b);
*/
/*for (i = 0; i < Nc; i++) {
	dK[i] = (phiL_f[i] / phiG_f[i] - phiL_b[i] / phiG_b[i]) / (2 * 0.001 * T);
}*/

		f = 0.0;
		df = 0.0;
		inner_counter = 0;
		for (i = 0; i < Nc; i++) {
			inner_counter++;
			f += ZI[i] * (K[i] - 1.0);
			//df += ZI[i] * dK[i];
			df += ZI[i] * K[i] * ((phiL_f[i] - phiL_b[i]) / (2 * 0.001 * P) - (phiG_f[i] - phiG_b[i]) / (2 * 0.001 * P));
		}
		P_old = P;
		P -= f / df;
		/*T_old = T;
		T -= f / df;*/
		//errT = abs((T - T_old) / T_old);

		for (i = 0; i < Nc; i++) {
			x[i] = ZI[i] / (1.0 + (K[i] - 1.0) * nv);
			y[i] = K[i] * x[i];
		}
		normalizArr(x, Nc);
		normalizArr(y, Nc);
		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiL);

		zeroInitializatinoCubicEOS(a, b, S, A, B, Z, phase);
		cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Z, phase);
		cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phiG);

		for (i = 0; i < Nc; i++) {
			LnphiG(i) = log(phiG[i]);
			LnphiL(i) = log(phiL[i]);
		}



		//auto tempvec = (LnAlpha0 - LnAlpha1).array() / LnAlpha0.array();
		//err_1 = tempvec.array().square().sum();
		//err_1 = tempvec.no
		//err_1 = abs((T_old - T) / T_old);
		err_1 = abs((P_old - P) / P_old);
	}

	return 0;
}




