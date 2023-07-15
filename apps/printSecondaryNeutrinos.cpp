#include "simprop.h"

using namespace simprop;

// void testProductionSpectrum() {
//   KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
//   utils::OutputFile out("test_neutrino_production_spectrum.txt");
//   out << "# x - spectrum\n";
//   out << std::scientific;
//   const auto units = SI::cm3 / SI::sec;
//   const auto epsCmb = 6.3e-4 * SI::eV;
//   const auto epsIr = 1e-2 * SI::eV;
//   auto EpAxis = utils::LogAxis<double>(1e18 * SI::eV, 1e22 * SI::eV, 201);
//   auto EnuAxis = utils::LogAxis<double>(1e15 * SI::eV, 1e22 * SI::eV, 201);
//   for (auto Enu : EnuAxis) {
//     for (auto Ep : EpAxis) {
//       auto x = Enu / Ep;
//       {
//         auto eta = 4. * epsCmb * Ep / pow2(SI::protonMassC2);
//         out << std::scientific << x << " " << eta << " ";
//         out << x * x * nuSpec.Phi(eta, x) / units << " ";
//       }
//       {
//         auto eta = 4. * epsIr * Ep / pow2(SI::protonMassC2);
//         out << x * x * nuSpec.Phi(eta, x) / units << " ";
//       }
//       out << "\n";
//     }
//   }
// }

double energydensity(std::shared_ptr<photonfields::CMB> cmb, size_t N = 8) {
  const auto epsMin = cmb->getMinPhotonEnergy();
  const auto epsMax = cmb->getMaxPhotonEnergy();
  auto integrand = [cmb](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    auto value = cmb->density(epsilon);
    return epsilon * value;
  };

  auto a = std::log(epsMin);
  auto b = std::log(epsMax);
  return utils::RombergIntegration<double>(integrand, a, b, N, 1e-4);
}

double I_deps(double Enu, double Ep, std::shared_ptr<photonfields::CMB> cmb,
              std::shared_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> nuSpec,
              size_t N = 8) {
  const auto epsMin = cmb->getMinPhotonEnergy();
  const auto epsMax = cmb->getMaxPhotonEnergy();

  auto integrand = [cmb, nuSpec, Enu, Ep](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = cmb->density(epsilon) * nuSpec->Phi(eta, Enu / Ep);
    return epsilon * value;
  };
  auto a = std::log(epsMin);
  auto b = std::log(epsMax);
  auto I = utils::RombergIntegration<double>(integrand, a, b, N, 1e-4);

  return I;
}

double testSourceTerm(double Enu, size_t N = 8) {
  const auto I_0 = 1.8e-30 * (1. / SI::eV / SI::m2 / SI::steradian / SI::sec);

  auto cmb = std::make_shared<photonfields::CMB>();
  auto nuSpec = std::make_shared<KelnerAharonian2008::NeutrinoProductionSpectrum>();

  auto integrand = [cmb, nuSpec, Enu](double lnEp) {
    auto Ep = std::exp(lnEp);
    const auto E_0 = 1e18 * SI::eV;
    return std::pow(Ep / E_0, -3.) * I_deps(Enu, Ep, cmb, nuSpec);
  };
  const auto lnEpMin = std::log(Enu);
  const auto lnEpMax = std::log(1e4 * Enu);
  auto I = utils::RombergIntegration<double>(integrand, lnEpMin, lnEpMax, N, 1e-3);

  return 4. * M_PI / SI::cLight * I_0 * I;
}

double testSourceTermDeltaFunction(double Enu, std::shared_ptr<photonfields::CMB> cmb,
                                   size_t N = 8) {
  const auto I_0 = 1.8e-30 * (1. / SI::eV / SI::m2 / SI::steradian / SI::sec);
  const auto epsilon = 6.3e-4 * SI::eV;

  const auto n_0 = energydensity(cmb);  //. / SI::cm3;
  LOGD << n_0 / (1. / SI::cm3);

  auto nuSpec = std::make_shared<KelnerAharonian2008::NeutrinoProductionSpectrum>();

  auto integrand = [epsilon, nuSpec, Enu](double lnEp) {
    auto Ep = std::exp(lnEp);
    const auto E_0 = 1e18 * SI::eV;
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    return std::pow(Ep / E_0, -3.) * nuSpec->Phi(eta, Enu / Ep);
  };
  const auto lnEpMin = std::log(Enu);
  const auto lnEpMax = std::log(1e4 * Enu);
  auto I = utils::RombergIntegration<double>(integrand, lnEpMin, lnEpMax, N, 1e-3);

  return 4. * M_PI / SI::cLight * I_0 * n_0 * I;
}

void testProductionSpectrum() {
  auto units = 1. / SI::eV / SI::m3 / SI::sec;
  auto EnuAxis = utils::LogAxis<double>(1e16 * SI::eV, 1e22 * SI::eV, 201);
  auto cmb = std::make_shared<photonfields::CMB>();

  utils::OutputFile out("test_neutrino_production_spectrum.txt");
  for (const auto& E_nu : EnuAxis) {
    out << std::scientific << E_nu / SI::eV << "\t";
    out << testSourceTerm(E_nu) / units << "\t";
    out << testSourceTermDeltaFunction(E_nu, cmb) / units << "\t";
    out << "\n";
  }
}

void testNeutrinoSpectrum(double gamma, double m, std::shared_ptr<photonfields::PhotonField> ebl,
                          std::string filename) {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto cmb = std::make_shared<photonfields::CMB>();

  std::vector<std::shared_ptr<losses::ContinuousLosses>> losses{
      std::make_shared<losses::PairProductionLosses>(cmb),
      std::make_shared<losses::PhotoPionContinuousLosses>(cmb)};
  auto b = solutions::Beniamino({gamma, m, 5e19 * SI::eV, 10.}, cosmology, losses).doCaching();
  auto nus_cmb = solutions::CosmoNeutrinos(b, cmb);
  auto nus_ebl = solutions::CosmoNeutrinos(b, ebl);
  {
    utils::OutputFile out(filename);
    const double units = SI::GeV / SI::cm2 / SI::sr / SI::sec;
    auto E = utils::LogAxis(1e14 * SI::eV, 1e20 * SI::eV, 8 * 5);
    for (const auto& E_nu : E) {
      out << std::scientific << E_nu / SI::eV << "\t";
      out << pow2(E_nu) * nus_cmb.computeNeutrinoFlux(E_nu, 6.) / units << "\t";
      out << pow2(E_nu) * nus_ebl.computeNeutrinoFlux(E_nu, 6.) / units << "\t";
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");

    // testProductionSpectrum();
    {
      auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
      testNeutrinoSpectrum(2.2, 0., ebl, "Saldana2021_2.2_0.txt");
      testNeutrinoSpectrum(2.4, 0., ebl, "Saldana2021_2.4_0.txt");
      testNeutrinoSpectrum(2.6, 0., ebl, "Saldana2021_2.6_0.txt");
      testNeutrinoSpectrum(2.4, -3., ebl, "Saldana2021_2.4_-3.txt");
      testNeutrinoSpectrum(2.4, 3., ebl, "Saldana2021_2.4_3.txt");
    }
    {
      auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
      testNeutrinoSpectrum(2.4, 0., ebl, "Dominguez2011_2.4_0.txt");
    }

    {
      auto ebl = std::make_shared<photonfields::Gilmore2012PhotonField>();
      testNeutrinoSpectrum(2.4, 0., ebl, "Gilmore2012_2.4_0.txt");
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}