#include "stdafx.h"
#include "CppUnitTest.h"

#include "global.h"

#define tol_unitConversion_tem 1e-6

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace Benchmark
{		
	TEST_CLASS(unit_conversion_global)
	{
	public:
		
		TEST_METHOD(F_to_C_test)
		{
			double F = 100.0;
			double Expected = 37.7777778;
			double Actual = F_to_C(F);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(F_to_K_test)
		{
			double F = 100.0;
			double Expected = 310.9277778;
			double Actual = F_to_K(F);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(F_to_R_test)
		{
			double F = 100.0;
			double Expected = 559.67;
			double Actual = F_to_R(F);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(C_to_F_test)
		{
			double C = 120.0;
			double Expected = 248.0;
			double Actual = C_to_F(C);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(C_to_K_test)
		{
			double C = 120.0;
			double Expected = 393.15;
			double Actual = C_to_K(C);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(C_to_R_test)
		{
			double C = 120.0;
			double Expected = 707.67;
			double Actual = C_to_R(C);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(R_to_C_test)
		{
			double R = 120.0;
			double Expected = -206.4833333;
			double Actual = R_to_C(R);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(R_to_K_test)
		{
			double R = 120.0;
			double Expected = 66.6666667;
			double Actual = R_to_K(R);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(R_to_F_test)
		{
			double R = 120.0;
			double Expected = -339.67;
			double Actual = R_to_F(R);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(K_to_C_test)
		{
			double K = 120.0;
			double Expected = -153.15;
			double Actual = K_to_C(K);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(K_to_F_test)
		{
			double K = 120.0;
			double Expected = -243.67;
			double Actual = K_to_F(K);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(K_to_R_test)
		{
			double K = 120.0;
			double Expected = 216.0;
			double Actual = K_to_R(K);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(m_to_foot_test)
		{
			double m = 120.0;
			double Expected = 393.7007874;
			double Actual = m_to_foot(m);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(m_to_inch_test)
		{
			double m = 120.0;
			double Expected = 4724.4094488;
			double Actual = m_to_inch(m);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(inch_to_meter_test)
		{
			double inch = 120.0;
			double Expected = 3.048;
			double Actual = inch_to_meter(inch);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(inch_to_foot_test)
		{
			double inch = 120.0;
			double Expected = 10.0;
			double Actual = inch_to_foot(inch);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(foot_to_meter_test)
		{
			double foot = 120.0;
			double Expected = 36.576;
			double Actual = foot_to_meter(foot);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(foot_to_inch_test)
		{
			double foot = 120.0;
			double Expected = 1440.0;
			double Actual = foot_to_inch(foot);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(bar_to_atm_test)
		{
			double bar = 120.0;
			double Expected = 118.43078032;
			double Actual = bar_to_atm(bar);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(bar_to_Pas_test)
		{
			double bar = 120.0;
			double Expected = 12000000.0;
			double Actual = bar_to_Pas(bar);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(bar_to_psi_test)
		{
			double bar = 120.0;
			double Expected = 1740.4528561;
			double Actual = bar_to_psi(bar);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(atm_to_bar_test)
		{
			double atm = 120.0;
			double Expected = 121.590012;
			double Actual = atm_to_bar(atm);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(atm_to_Pas_test)
		{
			double atm = 120.0;
			double Expected = 12159001.2;
			double Actual = atm_to_Pas(atm);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(atm_to_psi_test)
		{
			double atm = 120.0;
			double Expected = 1763.5140305;
			double Actual = atm_to_psi(atm);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(psi_to_atm_test)
		{
			double psi = 120.0;
			double Expected = 8.1655148477;
			double Actual = psi_to_atm(psi);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(psi_to_bar_test)
		{
			double psi = 120.0;
			double Expected = 8.273708736;
			double Actual = psi_to_bar(psi);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(psi_to_Pas_test)
		{
			double psi = 120.0;
			double Expected = 827370.8736;
			double Actual = psi_to_Pas(psi);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(Day_to_s_test)
		{
			double Day = 120.0;
			double Expected = 10368000;
			double Actual = Day_to_s(Day);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(s_to_Day_test)
		{
			double s = 120.0;
			double Expected = 0.0013888889;
			double Actual = s_to_Day(s);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scf_to_scm_test)
		{
			double scf = 12.0;
			double Expected = 0.339802159104;
			double Actual = scf_to_scm(scf);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scf_to_STB_test)
		{
			double scf = 12.0;
			double Expected = 2.13729128014842;
			double Actual = scf_to_STB(scf);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scm_to_scf_test)
		{
			double scm = 12.0;
			double Expected = 423.776000657863;
			double Actual = scm_to_scf(scm);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scm_to_gallon_test)
		{
			double scm = 12.0;
			double Expected = 3170.06462829778;
			double Actual = scm_to_gallon(scm);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(STB_to_scf_test)
		{
			double STB = 12.0;
			double Expected = 67.375;
			double Actual = STB_to_scf(STB);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(STB_to_scm_test)
		{
			double STB = 12.0;
			double Expected = 1.907847539136;
			double Actual = STB_to_scm(STB);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(STB_to_gallon_test)
		{
			double STB = 12.0;
			double Expected = 504.0;
			double Actual = STB_to_gallon(STB);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scfPerSTB_to_scmPerscm_test)
		{
			double scfPerSTB = 12.0;
			double Expected = 2.13729128014842;
			double Actual = scfPerSTB_to_scmPerscm(scfPerSTB);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(scmPerscm_To_scfPerSTB_test)
		{
			double scfPerSTB = 12.0;
			double Expected = 67.375;
			double Actual = scmPerscm_To_scfPerSTB(scfPerSTB);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(PaS_to_cp_test)
		{
			double PaS = 10.0;
			double Expected = 10000.0;
			double Actual = PaS_to_cp(PaS);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}

		TEST_METHOD(cp_to_PaS_test)
		{
			double cp = 10.0;
			double Expected = 0.01;
			double Actual = cp_to_PaS(cp);
			Assert::IsTrue(AlmostEqual(Expected, Actual, tol_unitConversion_tem));
			// TODO: Your test code here
		}


	};
}