#include "mathFunc.h"




void mathFunc::cubicEqSolver(double a1, double a2, double a3, vector<double>& ans) {
	double Q = (3 * a2 - a1 * a1) / 9.0;
	double J = (9.0 * a1 * a2 - 27.0 * a3 - 2 * a1 * a1 * a1) / 54.0;
	double D = Q * Q * Q + J * J;
	int rootNum;
	if (D > 0) {
		rootNum = 1;
		double pp = J - sqrt(D);
		ans.push_back(cbrt(J + sqrt(D)) + cbrt(J - sqrt(D)) - a1 / 3.0);	
	}
	else if (D < 0) {
		rootNum = 3;
		double theta = acos(J / sqrt(-Q * Q * Q));
		ans.push_back(2 * sqrt(-Q) * cos(theta / 3) - a1 / 3.0);
		ans.push_back(2 * sqrt(-Q) * cos(theta / 3.0 + 120.0 / 180.0 * piNum) - a1 / 3.0);
		ans.push_back(2 * sqrt(-Q) * cos(theta / 3.0 + 240.0 / 180.0 * piNum) - a1 / 3.0);
	}
	else {
		rootNum = 2;
		ans.push_back(2.0 * cbrt(J) - a1 / 3.0);
		ans.push_back (-cbrt(J) - a1 / 3.0);
	}

}

void mathFunc::matrixProduct(double** A, double** B, double** ANS, int N) {

	register int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			ANS[i][j] = vectorDotProduct(A[i], B[j], N);
		}
	}
}

double mathFunc::vectorDotProduct(double* a, double* b, int N) {
	double sum = 0;
	register int i;
	for (i = 0; i < N; i++)
		sum += a[i] * b[i];

	return sum;
}


void mathFunc::vectorA_CrossProduct_VectorAT(double* a, double* b, double** ANS, int N) {

	register int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			ANS[i][j] = a[i] * b[j];
		}
	}

}

// Calculation of cross product of matrix A(N,N) and vector b(N,1)
void mathFunc::matrixA_product_vectorb(double** A, double* b, double* ans, int N) {
	register int k, n;

	for (k = 0; k < N; k++) {
		ans[k] = 0.0;
	}

	for (k = 0; k < N; k++) {
		for (n = 0; n < N; n++) {
			ans[k] += A[k][n] * b[n];
		}
	}
}

void mathFunc::vectorb_product_matrixA(double** A, double* b, double* ans, int N) {
	register int k, n;
	for (k = 0; k < N; k++) {
			ans[k] = 0;
	}

	for (k = 0; k < N; k++) {
		for (n = 0; n < N; n++) {
			ans[k] += b[n] * A[n][k];
		}
	}
}

void mathFunc::matrixDividedByNumber(double** A, double c, int N) {
	register int k, n;

	for (k = 0; k < N; k++) {
		for (n = 0; n < N; n++) {
			A[k][n] /= c;
		}
	}
}


// Calculation of summation of two matrix A and B of size N
void mathFunc::matrixSummation(double** A, double** B, double** ANS, int N) {
	register int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			ANS[i][j] = A[i][j] + B[i][j];
		}
	}

}


void mathFunc::vectorSubtraction(double* a, double* b, double* ans, int N) {

	register int i;
	for (i = 0; i < N; i++) {
		ans[i] = a[i] - b[i];
	}
}

// Calculation of summation of two vector a and b of size N
void mathFunc::vectorSummation(double* a, double* b, double* ans, int N) {

	register int i;

	for (i = 0; i < N; i++) {
		ans[i] = a[i] + b[i];
	}

}


// Calculation of summation of  vector a of size N and constant c
void mathFunc::vectorSumWithNumber(double* a, double c, double* ans, int N) {

	register int i;

	for (i = 0; i < N; i++) {
		ans[i] = a[i] + c;
	}
}


// Calculation of product of  vector a of size N and constant c
void mathFunc::vectorProductbyNumber(double* a, double c, double* ans, int N) {

	register int i;

	for (i = 0; i < N; i++) {
		ans[i] = a[i] * c;
	}

}

double mathFunc::norm2Vector(double* a, int N) {

	register int i;
	double norm = 0;
	for (i = 0; i < N; i++) {
		norm += a[i] * a[i];
	}

	return norm;
}


// Calculation of error for vector x of size N
double mathFunc::errorSolution(double* xNew, double* xOld, int N) {

	register int i;
	double err = 0;
	for (i = 0; i < N; i++) {
		err += abs((xNew[i] - xOld[i]) / xOld[i]);
	}

	return err;
}



// Transpose of matrix A of size N
void mathFunc::matrixTranspose(double** A, double** AT, int N) {

	register int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			AT[i][j] = A[j][i];
		}
	}
}


void mathFunc::assignVecToVec(double* a1, double* a2, int N) {

	register int i;
	for (i = 0; i < N; i++) {
		a2[i] = a1[i];
	}
}


void mathFunc::assignMatToMat(double** A1, double** A2, int N) {

	register int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A2[i][j] = A1[i][j];
		}
	}


}





mathFunc::mathFunc()
{
}


mathFunc::~mathFunc()
{
}
