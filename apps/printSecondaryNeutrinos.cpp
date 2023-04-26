#include "simprop.h"

using namespace simprop;

void testProductionSpectrum() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  utils::OutputFile out("test_neutrino_production_spectrum.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  const auto epsCmb = 6.3e-4 * SI::eV;
  const auto epsIr = 1e-2 * SI::eV;
  auto EpAxis = utils::LogAxis<double>(1e18 * SI::eV, 1e22 * SI::eV, 201);
  auto EnuAxis = utils::LogAxis<double>(1e15 * SI::eV, 1e22 * SI::eV, 201);
  for (auto Enu : EnuAxis) {
    for (auto Ep : EpAxis) {
      auto x = Enu / Ep;
      {
        auto eta = 4. * epsCmb * Ep / pow2(SI::protonMassC2);
        out << std::scientific << x << " " << eta << " ";
        out << x * x * nuSpec.Phi(eta, x) / units << " ";
      }
      {
        auto eta = 4. * epsIr * Ep / pow2(SI::protonMassC2);
        out << x * x * nuSpec.Phi(eta, x) / units << " ";
      }
      out << "\n";
    }
  }
}

void testNeutrinoSpectrum() {
  solutions::CosmoNeutrinos nu(2.6, 0., -1.);
  {
    utils::OutputFile out("SimProp_proton_spectrum.txt");
    const double units = pow2(SI::eV) / SI::m2 / SI::sr / SI::sec;
    auto E = utils::LogAxis(1e16 * SI::eV, 1e21 * SI::eV, 100);
    for (const auto& E_p : E) {
      out << std::scientific << E_p / SI::eV << "\t";
      out << pow3(E_p) * nu.getProtonFlux(E_p, 0.) / units << "\t";
      out << pow3(E_p) * nu.getProtonFlux(E_p, 3.) / units << "\t";
      out << pow3(E_p) * nu.getProtonFlux(E_p, 6.) / units << "\t";
      out << "\n";
    }
  }
  // {
  //   utils::OutputFile out("SimProp_neutrino_spectrum.txt");
  //   const double units = SI::GeV / SI::cm2 / SI::sr / SI::sec;
  //   auto E = utils::LogAxis(1e16 * SI::eV, 1e20 * SI::eV, 8 * 4);
  //   for (const auto& E_nu : E) {
  //     out << std::scientific << E_nu / SI::eV << "\t";
  //     out << E_nu * E_nu * nu.computeNeutrinoFlux(E_nu, 5.) / units << "\t";
  //     out << "\n";
  //   }
  // }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    // testProductionSpectrum();
    testNeutrinoSpectrum();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}