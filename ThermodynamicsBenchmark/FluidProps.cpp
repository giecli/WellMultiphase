#include "stdafx.h"
#include "CppUnitTest.h"



using namespace Microsoft::VisualStudio::CppUnitTestFramework;

#include"global.h"
#include"BlackOil.h"
#include"Compositional.h"



using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace ThermodynamicsBenchmark
{		
	TEST_CLASS(BlackOilTest)
	{
	public:
		
		TEST_METHOD(Pb_Standing)
		{
			BlackOil sample;
			sample.set_T(F_to_K(220));
			sample.set_gammaGas(0.768);
			sample.set_API(40.7);
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			double Expected = 2738.4140311201;
			sample.bubblePointPressure("Standing");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Pb_Lasater)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(120));
			sample.set_gammaOil(0.8353);
			sample.set_gammaGas(0.804);
			sample.set_API(37.9);
			sample.set_MwOil(193.0);
			sample.set_Rs(scfPerSTB_to_scmPerscm(769.0));
			double Expected = 2050;
			sample.bubblePointPressure("Lasater");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Pb_Glaso)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(120));
			sample.set_gammaGas(0.804);
			sample.set_API(40.7);
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			double Expected = 2647.56287968155;
			sample.bubblePointPressure("Glaso");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		
		TEST_METHOD(Pb_Beggs_Vazquez1)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(120));
			sample.set_gammaGas(0.804);
			sample.set_API(25.0);
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			double Expected = 4000.05663050777;
			sample.bubblePointPressure("Beggs-Vazquez");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Pb_Beggs_Vazquez2)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(120));
			sample.set_gammaGas(0.804);
			sample.set_API(40.7);
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			double Expected = 2343.18759589129;
			sample.bubblePointPressure("Beggs-Vazquez");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Pb_Dindoruk_Christman)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(120));
			sample.set_gammaGas(0.804);
			sample.set_API(35.7);
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			double Expected = 3133.68575628165;
			sample.bubblePointPressure("Dindoruk-Christman");
			double Actual = Pas_to_psi(sample.get_Pb());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Standing)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_gammaGas(0.786);
			sample.set_API(40.7);
			sample.set_Pb(psi_to_Pas(2685.0));
			double Expected = 768.0;
			sample.SaturatedGasOilRatio("Standing");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Lasater)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_gammaGas(0.786);
			sample.set_API(40.7);
			sample.stockTankMw("PVTsim");
			sample.set_gammaOil(0.90);
			sample.set_Pb(psi_to_Pas(2685.0));
			double Expected = 768.0;
			sample.SaturatedGasOilRatio("Lasater");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Vazquez_Beggs1)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_gammaGas(0.786);
			sample.set_Pb(psi_to_Pas(2685.0));
			sample.set_API(25.0);
			double Expected = 689.221990429723;
			sample.SaturatedGasOilRatio("Vazquez-Beggs");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Vazquez_Beggs2)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_gammaGas(0.786);
			sample.set_Pb(psi_to_Pas(2685.0));
			sample.set_API(40.7);
			double Expected = 412.377087147648;
			sample.SaturatedGasOilRatio("Vazquez-Beggs");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Glaso)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Pb(psi_to_Pas(2685.0));
			sample.set_gammaGas(0.786);
			sample.set_API(40.7);
			double Expected = 672.160499847815;
			sample.SaturatedGasOilRatio("Glaso");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Rs_Dindoruk_Christman)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Pb(psi_to_Pas(2685.0));
			sample.set_gammaGas(0.786);
			sample.set_API(25.0);
			double Expected = 465.074969879656;
			sample.SaturatedGasOilRatio("Dindoruk-Christman");
			double Actual = scmPerscm_To_scfPerSTB(sample.get_Rs());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1));
			// TODO: Your test code here
		}

		TEST_METHOD(Bob_Standing)
		{
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_gammaOil(0.904);
			
			double Expected = 1.45923277766733;
			sample.OFVF("Standing");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(Bob_Glaso)
		{
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_gammaOil(0.904);

			double Expected = 1.41823342358221;
			sample.OFVF("Glaso");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(Bob_AlMarhoun)
		{
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_gammaOil(0.904);

			double Expected = 1.41505306488191;
			sample.OFVF("AlMarhoun");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(Bob_Vazquez_Beggs1)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_API(25.0);

			double Expected = 1.51908408346056;
			sample.OFVF("Vazquez-Beggs");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}


		TEST_METHOD(Bob_Vazquez_Beggs2)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_API(35.0);

			double Expected = 1.51908408346056;
			sample.OFVF("Vazquez-Beggs");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(Bob_Dindoruk_Christman)
		{
			
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_Rs(scfPerSTB_to_scmPerscm(768.0));
			sample.set_gammaGas(0.786);
			sample.set_API(35.0);
			sample.gammaOilCal();

			double Expected = 1.42140076492475;
			sample.OFVF("Dindoruk-Christman");
			double Actual = sample.get_Bob();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDead_Beal_Standing)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_API(35.0);

			double Expected = 1000.0;
			sample.DeadOilViscosity("Beal-Standing");
			double Actual = PaS_to_cp(sample.get_muDeadOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDead_Beggs_Robinson)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_API(35.0);

			double Expected = 1000.0;
			sample.DeadOilViscosity("Beggs-Robinson");
			double Actual = PaS_to_cp(sample.get_muDeadOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDead_Glaso)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_API(35.0);

			double Expected = 1000.0;
			sample.DeadOilViscosity("Glaso");
			double Actual = PaS_to_cp(sample.get_muDeadOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDead_AlKhafaji)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_API(35.0);

			double Expected = 1000.0;
			sample.DeadOilViscosity("AlKhafaji");
			double Actual = PaS_to_cp(sample.get_muDeadOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDead_Dindoruk_Christman)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(220.0));
			sample.set_API(35.0);
			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_Pb(psi_to_Pas(2000.0));

			double Expected = 1000.0;
			sample.DeadOilViscosity("Dindoruk-Christman");
			double Actual = PaS_to_cp(sample.get_muDeadOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muSat_Standing)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 1000.0;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muSat_Beggs_Robinson)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 0.342263060707993;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muSat_Bergman)
		{														//********************************
			// Could not find any test in the literature, but the answer *****DOES NOT seem**** OK
			BlackOil sample;									//********************************

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 1000.0;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muSat_Aziz)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 0.342263060707993;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muSat_AlKhafaji)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 0.342263060707993;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(muDindoruk_Christman)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;

			sample.set_Rs(scfPerSTB_to_scmPerscm(786.0));
			sample.set_muDeadOil(cp_to_PaS(1.1));

			double Expected = 0.342263060707993;
			sample.SatOilViscosity();
			double Actual = PaS_to_cp(sample.get_muSatOil());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 0.1));
			// TODO: Your test code here
		}

		TEST_METHOD(Gas_FVF)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_P(psi_to_Pas(2000.0));
			sample.set_T(R_to_K(520.0));
			sample.set_Z(1.2);
			double Expected = 0.0088231416;
			double Actual = sample.BgCalculation();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}
		
		TEST_METHOD(Oil_gammaOil)
		{
			//
			BlackOil sample;
			sample.set_API(12);
			
			double Expected = 0.9860627177700348;
			double Actual = sample.gammaOilCal();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Gas_Pc_Tc_Sutton)
		{
			
			BlackOil sample;
			sample.set_gammaGas(0.9);

			double ExpectedPpc = 635.984;
			double ExpectedTpc = 423.81;
			sample.gasCriticalProps("Sutton");
			double ActualPpc = Pas_to_psi(sample.get_gasPpc());
			double ActualTpc = K_to_R(sample.get_gasTpc());
			Assert::IsTrue(AlmostEqual(ExpectedPpc, ActualPpc, 1e-6));
			Assert::IsTrue(AlmostEqual(ExpectedTpc, ActualTpc, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Oil_stockTankMw)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_API(12);

			double Expected = 0.9860627177700348;
			double Actual = sample.stockTankMw("PVTsim");
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Gas_viscosity_Dempsey)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(250.0));
			sample.set_gammaGas(0.9);
			sample.set_P(psi_to_Pas(1200.0));
			sample.gasMolecularWeight();
			double Expected = 0.9860627177700348;
			sample.GasViscosity("Dempsey");
			double Actual = PaS_to_cp(sample.get_muGas());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Gas_viscodity_Lee_Gunzalez)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(250.0));
			sample.set_gammaGas(0.9);
			sample.set_P(psi_to_Pas(1200.0));
			sample.gasMolecularWeight();
			double Expected = 0.9860627177700348;
			sample.GasViscosity("Lee-Gunzalez");
			double Actual = PaS_to_cp(sample.get_muGas());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}


		TEST_METHOD(Gas_Lucas)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(250.0));
			sample.set_gammaGas(0.9);
			sample.set_P(psi_to_Pas(1200.0));
			sample.gasCriticalProps("Sutton");
			sample.gasMolecularWeight();
			double Expected = 0.9860627177700348;
			//sample.ZFactor("PVTsim");
			sample.GasViscosity("Lucas");
			double Actual = PaS_to_cp(sample.get_muGas());
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

		TEST_METHOD(Gas_Zfactor)
		{
			// Could not find any test in the literature, but the answer seems OK
			BlackOil sample;
			sample.set_T(F_to_K(250.0));
			sample.set_gammaGas(0.9);
			sample.set_P(psi_to_Pas(6900.0));
			sample.gasCriticalProps("Sutton");
			double Expected = 0.9860627177700348;
			sample.ZFactor();
			double Actual = sample.get_Z();
			Assert::IsTrue(AlmostEqual(Expected, Actual, 1e-6));
			// TODO: Your test code here
		}

	};

	TEST_CLASS(CompositionalModel)
	{
	public:

		TEST_METHOD(RK_TEST)
		{
			double T = 300.0;
			double P = atm_to_Pas(11.81);
			Compositional sample;
			int N = 1;
			sample.setNc(N);
			sample.setEOS("RK");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}
			
			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				BIC[i][i] = 0;
			}
			Tc[0] = 365.0;
			Pc[0] = atm_to_Pas(46.5);
			comp[0] = 1.0;
			sample.RKEOS(T, P, Tc, Pc, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 0.8255;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(SRK_TEST_mixingRules)
		{
			double T = C_to_K(70.0);
			double P = bar_to_Pas(200.0);
			Compositional sample;
			int N = 2;
			sample.setNc(N);
			sample.setEOS("SRK");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}
			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			comp[0] = 0.4;
			comp[1] = 0.6;
			ACF[1] = 0.09800;
			sample.SRKEOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 0.7319;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_mixingRules)
		{
			double T = C_to_K(70.0);
			double P = bar_to_Pas(200.0);
			Compositional sample;
			int N = 2;
			sample.setNc(N);
			sample.setEOS("PR");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}
			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			comp[0] = 0.4;
			comp[1] = 0.6;
			ACF[1] = 0.09800;
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 0.6781;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR78_TEST_mixingRules)
		{
			double T = C_to_K(70.0);
			double P = bar_to_Pas(200.0);
			Compositional sample;
			int N = 2;
			sample.setNc(N);
			sample.setEOS("PR78");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}
			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			comp[0] = 0.4;
			comp[1] = 0.6;
			ACF[1] = 0.09800;
			sample.PR78EOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 0.6781;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR78_TEST_mixingRules2)
		{
			double T = C_to_K(70.0);
			double P = bar_to_Pas(200.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			sample.setEOS("PR78");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;
			sample.PR78EOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 1.3210;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_fugacityCoef)
		{
			double T = 477.6;
			double P = bar_to_Pas(18.6);
			Compositional sample;
			int N = 1;
			sample.setNc(N);
			sample.setEOS("PR");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* phi = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				phi[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 507.6;
			Pc[0] = bar_to_Pas(30.25);
			ACF[0] = 0.3013;
			comp[0] = 1.0;
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			sample.cubicFugacityCoef(A, B, a, b, bi, ai, S, Z, phi);
			double Exoected = 0.729704;
			Assert::IsTrue(AlmostEqual(Exoected, phi[0], 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR78_TEST_fugacity)
		{
			double T = C_to_K(70.0);
			double P = bar_to_Pas(200.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			sample.setEOS("PR78");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;
			sample.PR78EOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 1.3210;
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_2SS)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double* x = new double[N];
			double* y = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				K[i] = 0;
				x[i] = 0;
				y[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;
			
			
			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.5;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			sample.twoPhaseFlash_SS(P, T, ai, bi, BIC, comp, K, nv, x, y);
			double Expected = 1.3210;
			Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_2MSS)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			sample.setEOS("PR78");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double* x = new double[N];
			double* y = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				K[i] = 0;
				x[i] = 0;
				y[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.5;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			sample.twoPhaseFlash_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
			double Exoected = 0.2010;
			Assert::IsTrue(AlmostEqual(Exoected, nv, 1e-4));
			// TODO: Your test code here
		}


		TEST_METHOD(PR_TEST_flashCalculation_QNN_MSS)
		{
			double T = C_to_K(200.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double* x = new double[N];
			double* y = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				K[i] = 0;
				x[i] = 0;
				y[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.8;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			sample.twoPhaseFlash_QNN_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
			double Exoected = 0.2272;
			Assert::IsTrue(AlmostEqual(Exoected, nv, 1e-4));
			// TODO: Your test code here
		}


		TEST_METHOD(PR_TEST_flashCalculation_QNN_Stability)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				K[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.8;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double Acctual = sample.twoPhaseStabilityFlash(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 1.0;
			Assert::IsTrue(Acctual >= Expected);
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_BFGS_Stability)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				K[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.8;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double Acctual = sample.twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 1.0;
			Assert::IsTrue(Acctual >= Expected);
			// TODO: Your test code here
		}



		TEST_METHOD(PR_TEST_flashCalculation_QNN_MSS_Stability)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.8;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double Actual = sample.twoPhaseStabilityFlash(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.2010;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if (Actual >= stat) {
				sample.twoPhaseFlash_QNN_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}

			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_BFGS_MSS_Stability)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.8;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double Actual = sample.twoPhaseStabilityFlash(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.2010;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if (Actual >= stat) {
				sample.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}

			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_BFGS_MSS_Stability_FullyImplicit)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.5;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double Actual = sample.twoPhaseStabilityFlash_BFGS_FullyImplicit(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.2025;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if (Actual >= stat) {
				nv = 0.2;
				//sample.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				sample.twoPhaseFlash_BFGS_FullyImplicit(P, T, ai, bi, BIC, comp, K, nv, x, y);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}

			// TODO: Your test code here
		}

		
		
		TEST_METHOD(Enthalpy_test)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(100.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			sample.setEOS("PR");
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double** BIC = new double*[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0;
				bi[i] = 0;
				S[i] = 0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.67011;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;
			sample.PR78EOS(T, P, Tc, Pc, ACF, ai, bi);
			string phase;
			double a, b, A, B, Z;
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			double Exoected = 1.3210;
			sample.setNc(1);
			double* phi = new double[1];
			sample.cubicFugacityCoef(A, B, ai[0], bi[0], bi, ai, ai, Z, phi);
			sample.setNc(3);
			double* a0 = new double[N];
			double* a1 = new double[N];
			double* a2 = new double[N];
			double* a3 = new double[N];
			double* a4 = new double[N];
			for (int i = 0; i < N; i++) {
				a0[i] = 0.0;
				a1[i] = 0.0;
				a2[i] = 0.0;
				a3[i] = 0.0;
				a4[i] = 0.0;
			}
			a0[0] = 4.568;
			a1[0] = -8.975;
			a2[0] = 3.631;
			a3[0] = -3.407;
			a4[0] = 1.091;
			
			a0[1] = 4.178;
			a1[1] = -4.427;
			a2[1] = 5.66;
			a3[1] = -6.651;
			a4[1] = 2.487;

			a0[2] = 21.18;
			a1[2] = -8.424;
			a2[2] = 39.969;
			a3[2] = 53.29;
			a4[2] = 21.482;

			double H = sample.enthalpyReal(comp, T, P, a0, a1, a2, a3, a4, Z, A, B, *phi);
			Assert::IsTrue(AlmostEqual(Exoected, Z, 1e-4));
			// TODO: Your test code here
		}

		TEST_METHOD(HeatCapacity_test)
		{
			double T = C_to_K(150.0);
			double P = bar_to_Pas(100.0);
			Compositional sample;
			int N = 2;
			sample.setNc(N);

			double* comp = new double[N];
			comp[0] = 0.8;
			comp[1] = 0.2;

			double* a0 = new double[N];
			double* a1 = new double[N];
			double* a2 = new double[N];
			double* a3 = new double[N];
			double* a4 = new double[N];
			for (int i = 0; i < N; i++) {
				a0[i] = 0.0;
				a1[i] = 0.0;
				a2[i] = 0.0;
				a3[i] = 0.0;
				a4[i] = 0.0;
			}
			a0[0] = 4.568;
			a1[0] = -8.975;
			a2[0] = 3.631;
			a3[0] = -3.407;
			a4[0] = 1.091;

			a0[1] = 4.178;
			a1[1] = -4.427;
			a2[1] = 5.66;
			a3[1] = -6.651;
			a4[1] = 2.487;


			double cp = sample.HeatCapacityIdeal(comp, T, a0, a1, a2, a3, a4);
			Assert::IsTrue(AlmostEqual(2, 3.0, 1e-4));
			// TODO: Your test code here
		}


		TEST_METHOD(LBC_TEST_vicosity)
		{
			double T = C_to_K(400.0);
			double P = bar_to_Pas(10.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* Vc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			double* Mw = new double[N];
			
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			Vc[0] = 0.09900 / 1000.0;
			ACF[0] = 0.00800;
			Mw[0] = 16.04 * 1e-3;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Mw[1] = 30.07 * 1e-3;
			Vc[1] = 0.14800 / 1000.0;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			Vc[2] = 0.85215 / 1000.0;
			ACF[2] = 0.67011;
			Mw[2] = 206.0 * 1e-3;

			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			
			sample.setEOS("PR");
			string phase;
			double Expected;
			double a, b, A, B, Z;
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			sample.cubicEOS(T, P, comp, ai, bi, a, b, S, A, B, BIC, Z, phase);
			Expected = 132.82 * 1e-3;
			double MwMix = sample.molecularWight(comp, Mw);
			Expected = 132.82 * 1e-3;
			Assert::IsTrue(AlmostEqual(Expected, MwMix, 1e-4));
			double d = sample.density_f(P, T, MwMix, Z);
			//Expected = 0.6416 * 1e+3;
			//Assert::IsTrue(AlmostEqual(Expected, d, 1.0));
			double eta = sample.viscosity_LBC(P, T, d, comp, Tc, Pc, Vc, Mw);
			Expected = 0.2681;
			Assert::IsTrue(AlmostEqual(Expected, eta, 1e-1));
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_flashCalculation_IFT)
		{
			double T = C_to_K(10.0);
			double P = bar_to_Pas(50.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			double* Par = new double[N];
			double* Mw = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Par[0] = 77.300;
			Mw[0] = 16.0429 * 1e-3;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Par[1] = 108.900;
			Mw[1] = 30.0698 * 1e-3;
			Tc[2] = 676.266;
			Pc[2] = bar_to_Pas(18.2409);
			ACF[2] = 0.7678;
			Par[2] = 541.340;
			Mw[2] = 206.00 * 1e-3;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			double a, b, A, B, Zliq, Zgas;
			double nv = 0.9;
			sample.setEOS("SRK");
			sample.SRKEOS(T, P, Tc, Pc, ACF, ai, bi);
			double Actual = sample.twoPhaseStabilityFlash_BFGS(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.0009078;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if ((Actual >= stat)) {
				sample.twoPhaseFlash_QNN_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}

			double MwGas = sample.molecularWight(y, Mw);
			double MwLiq = sample.molecularWight(x, Mw);
			Zliq = 0.0;
			Zgas = 0.0;
			A = 0.0;
			B = 0.0;
			a = 0.0;
			b = 0.0;
			//sample.zeroInitializatinoCubicEOS(a, b, S, A, B, Zliq, phase);
			sample.cubicEOS(T, P, x, ai, bi, a, b, S, A, B, BIC, Zliq, phase);
			sample.zeroInitializatinoCubicEOS(a, b, S, A, B, Zliq, phase);
			sample.cubicEOS(T, P, y, ai, bi, a, b, S, A, B, BIC, Zgas, phase);

			double rhoL = sample.density_f(P, T, MwLiq, Zliq);
			double rhoG = sample.density_f(P, T, MwGas, Zgas);
			double IFT = sample.GasOilIFT(rhoL / (MwLiq * 1e+6), rhoG / (MwGas * 1e+6), Par, x, y);
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_saturationT_BFGS_MSS_Stability_bubble)
		{
			double T = C_to_K(-40.46);
			double P = bar_to_Pas(20.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.0;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			int Tsat = sample.twoPhaseFlashEnv_BFGS_MSS_bubble(P, T, ai, bi, BIC, comp, nv, x, y, Tc, Pc, ACF);

			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_saturationT_BFGS_MSS_Stability_dew)
		{
			double T = C_to_K(280.0);
			double P = bar_to_Pas(5.0);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 1.0;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			int Tsat = sample.twoPhaseFlashEnv_BFGS_MSS_dew(P, T, ai, bi, BIC, comp, nv, x, y, Tc, Pc, ACF);

			// TODO: Your test code here
		}



		TEST_METHOD(PR_TEST_saturationT_BFGS_MSS_Stability_new)
		{
			double T = C_to_K(2.85000610351563);
			double P = bar_to_Pas(42.1215362548828);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.0;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double T_old = C_to_K(-40.46);
			double P_old0 = bar_to_Pas(23.13);
			
			double Actual = sample.twoPhaseStabilityFlash_BFGS_FullyImplicit(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.2025;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if (Actual >= stat) {
				nv = 0.2;
				//sample.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				sample.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}
			// TODO: Your test code here
		}

		TEST_METHOD(PR_TEST_saturationT_BFGS_MSS_Stability_dew_newTest)
		{
			double T = C_to_K(2.85000610351563);
			double P = bar_to_Pas(42.1215362548828);
			Compositional sample;
			int N = 3;
			sample.setNc(N);
			double* ai = new double[N];
			double* bi = new double[N];
			double* Tc = new double[N];
			double* Pc = new double[N];
			double* comp = new double[N];
			double* S = new double[N];
			double* ACF = new double[N];
			double* K = new double[N];
			double** BIC = new double*[N];
			double* x = new double[N];
			double* y = new double[N];
			for (int i = 0; i < N; i++) {
				BIC[i] = new double[N];
			}

			for (int i = 0; i < N; i++) {
				ai[i] = 0.0;
				bi[i] = 0.0;
				S[i] = 0.0;
				K[i] = 0.0;
				x[i] = 0.0;
				for (int j = 0; j < N; j++) {
					BIC[i][j] = 0;
				}
			}

			Tc[0] = 190.60;
			Pc[0] = bar_to_Pas(46.0015);
			ACF[0] = 0.00800;
			Tc[1] = 305.40;
			Pc[1] = bar_to_Pas(48.8387);
			ACF[1] = 0.09800;
			Tc[2] = 697.973;
			Pc[2] = bar_to_Pas(17.8578);
			ACF[2] = 0.6863;
			comp[0] = 0.2;
			comp[1] = 0.2;
			comp[2] = 0.6;


			sample.wilsonK(P, T, Tc, Pc, ACF, K);
			string phase;
			//double a, b, A, B, Z;
			double nv = 0.0;
			sample.setEOS("PR");
			sample.PREOS(T, P, Tc, Pc, ACF, ai, bi);
			double T_old = C_to_K(-40.46);
			double P_old0 = bar_to_Pas(23.13);

			double Actual = sample.twoPhaseStabilityFlash_BFGS_FullyImplicit(P, T, ai, bi, BIC, comp, K, nv);
			double Expected = 0.2025;
			Assert::IsTrue(Actual >= Expected);
			double stat = 1.0;
			if (Actual >= stat) {
				nv = 0.2;
				//sample.twoPhaseFlash_BFGS_MSS(P, T, ai, bi, BIC, comp, K, nv, x, y);
				sample.twoPhaseFlashEnv_BFGS_MSS_dew(P, T, ai, bi, BIC, comp, nv, x, y, Tc, Pc, ACF);
				Assert::IsTrue(AlmostEqual(Expected, nv, 1e-4));
			}
			// TODO: Your test code here
		}




	};
}
