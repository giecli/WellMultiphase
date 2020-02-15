#include "stdafx.h"
#include "CppUnitTest.h"
#include "global.h"
#include "Preliminary.h"
#include "twoPhaseModels.h"
#include "BlackOil.h"
#include "Compositional.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace PipeFlowTest
{
	TEST_CLASS(UnitTest1)
	{
	public:

		TEST_METHOD(TestMethod1)
		{
			double v_sL = foot_to_meter(1.601);
			double mu_L = cp_to_PaS(13.09);
			double v_sg = foot_to_meter(2.784);
			double mu_g = cp_to_PaS(0.019);
			double rho_L = bmft3_to_kgm3(55.042);
			double rho_g = bmft3_to_kgm3(2.17);
			double sigma = 31.6 * 1e-3;
			double d = inch_to_meter(2.99);
			twoPhaseModels ansari;
			Preliminary sample;
			double f_g = v_sg / (v_sg + v_sL);
			double v_m = v_sg + v_sL;
			double Re_sL = sample.ReynoldNumber(rho_L, v_sL, d, mu_L);
			double f_sL = sample.colebrookWhiteMoody(d * 1e-4, d, Re_sL);
			string pattern = ansari.ansariPattern(rho_L, rho_g, sigma, d, f_g, v_m, v_sg, f_sL);
			double f_L, f_TP, rho_m, mu_m, dv_mdz, e;
			e = d * 1e-4;
			ansari.ansariPressureGradient(pattern, f_L, f_TP, rho_m, mu_m, v_m, d, piNum / 2.0, rho_L, rho_g, mu_L, mu_g, v_sL,
				v_sg, dv_mdz, e, sigma);
			//TODO: Your test code here
		}

		TEST_METHOD(FLOW_PROFILE_oil_AND_gas) {
			double WellDepth = foot_to_meter(5151.0);
			double Qo = STB_to_scm(1140.0);
			double Pw = psi_to_Pas(505.0);
			int gridNum = 50;
			double dL = WellDepth / gridNum;
			double Tw = F_to_K(93.0);
			double P_boundary = psi_to_Pas(505.0);
			double gammaGas = 0.8;
			double API = 23.0;
			double Rs = scfPerSTB_to_scmPerscm(450.0);
			double Qg = Rs * Qo;
			double d = inch_to_meter(2.99);
			double e = d * 1e-4;
			double theta = 90.0 / 180.0 * piNum;
			double Co = (1e-5) / psi_to_Pas(1.0);
			string units = "SI";
			string Pseudocritical_method = "T2";
			string MwMethod = "PVTsim";
			string muGasMethod = "Lee-Gunzalez";
			string muDeadOilMethod = "Beggs-Robinson";
			string muSatOilMethod = "Beggs-Robinson";
			string ZFactorMethod = "PVTsim";
			string BobMethod = "Standing";
			string flowDir = "upstream";
			twoPhaseModels FlowProfileObj;
			double* PAns = new double[gridNum + 1];
			double* L = new double[gridNum + 1];
			FlowProfileObj.FlowProfile(Qg, Qo, Pw, gridNum, Tw, Pw, gammaGas, API,
				Rs, Co, e, d, dL, theta, "SI", Pseudocritical_method,
				MwMethod, muGasMethod, muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir,
				PAns, L);
		}

		TEST_METHOD(FLOW_PROFILE_oil_AND_gas_AND_water) {
			double WellDepth = foot_to_meter(5355.0);
			double Qo = STB_to_scm(59.0);
			double Qw = STB_to_scm(542.0);
			double waterSal = 35000 * 1e-6;
			double Pw = psi_to_Pas(113.0);
			int gridNum = 20;
			double dL = WellDepth / gridNum;
			double Tw = F_to_K(85.0);
			double P_boundary = psi_to_Pas(113.0);
			double gammaGas = 1.04;
			double API = 34.3;
			double Qg = scf_to_scm(41000.0);
			double Rs = Qg / Qo;
			double d = foot_to_meter(0.1198 * 2.0);
			double e = 2.0 * d * 1e-4;
			double theta = 90.0 / 180.0 * piNum;
			double Co = (1e-5) / psi_to_Pas(1.0);
			string units = "SI";
			string Pseudocritical_method = "T2";
			string MwMethod = "PVTsim";
			string muGasMethod = "Lee-Gunzalez";
			string muDeadOilMethod = "Beggs-Robinson";
			string muSatOilMethod = "Beggs-Robinson";
			string ZFactorMethod = "PVTsim";
			string BobMethod = "Standing";
			string flowDir = "upstream";
			twoPhaseModels FlowProfileObj;
			double* PAns = new double[gridNum + 1];
			double* L = new double[gridNum + 1];
			FlowProfileObj.FlowProfile(Qg, Qo, Qw, waterSal, Pw, gridNum, Tw, Pw, gammaGas, API,
				Rs, Co, e, d, dL, theta, "SI", Pseudocritical_method,
				MwMethod, muGasMethod, muDeadOilMethod, muSatOilMethod, ZFactorMethod, BobMethod, flowDir, PAns, L);
		}


	};
}