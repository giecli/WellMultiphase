#pragma once


#include<vector>
#include<cmath>
#include"global.h"


using namespace std;

class mathFunc
{
public:
	mathFunc();
	~mathFunc();
	void cubicEqSolver(double a1, double a2, double a3, vector<double>& ans);
	void matrixProduct(double** A, double** B, double** ANS, int N);
	double vectorDotProduct(double* a, double* b, int N);
	void vectorA_CrossProduct_VectorAT(double* a, double* b, double** ANS, int N);
	void matrixA_product_vectorb(double** A, double* b, double* ans, int N);
	void vectorb_product_matrixA(double** A, double* b, double* ans, int N);
	void matrixSummation(double** A, double** B, double** ANS, int N);
	void vectorSubtraction(double* a, double* b, double* ans, int N);
	void vectorSummation(double* a, double* b, double* ans, int N);
	void vectorSumWithNumber(double* a, double c, double* ans, int N);
	void vectorProductbyNumber(double* a, double c, double* ans, int N);
	double norm2Vector(double* a, int N);
	double errorSolution(double* xNew, double* xOld, int N);
	void matrixTranspose(double** A, double** AT, int N);
	void assignVecToVec(double* a1, double* a2, int N);
	void assignMatToMat(double** A1, double** A2, int N);
	void matrixDividedByNumber(double** A, double c, int N);
	


};

