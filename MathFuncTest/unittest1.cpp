#include "stdafx.h"
#include "CppUnitTest.h"
#include"mathFunc.h"
#include<vector>
using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MathFuncTest
{		
	TEST_CLASS(MathFuncTests)
	{
	public:
		
		TEST_METHOD(cubicEqSolver_test1)
		{
			mathFunc sample;
			double a1 = -1.90257;
			double a2 = +1.55625;
			double a3 = -1.40463;
			vector<double> ans;
			sample.cubicEqSolver(a1, a2, a3, ans);
			double ExpectedAns = 1.49069;
			int ExpectedRootNum = 1;
			Assert::IsTrue(AlmostEqual(ExpectedAns, ans[0], 1e-4));
			Assert::IsTrue(ExpectedRootNum == ans.size());

			// TODO: Your test code here
		}

		TEST_METHOD(cubicEqSolver_test2)
		{
			mathFunc sample;
			double a1 = -0.9491604;
			double a2 = +0.22381034;
			double a3 = -0.0142259;
			vector<double> ans;
			sample.cubicEqSolver(a1, a2, a3, ans);
			vector<double> ExpectedAns = { 0.62954, 0.10557, 0.21405 };
			Assert::IsTrue(AlmostEqual(ExpectedAns[0], ans[0], 1e-4));
			Assert::IsTrue(AlmostEqual(ExpectedAns[1], ans[1], 1e-4));
			Assert::IsTrue(AlmostEqual(ExpectedAns[2], ans[2], 1e-4));
			Assert::IsTrue(ExpectedAns.size() == ans.size());

			// TODO: Your test code here
		}


	};
}