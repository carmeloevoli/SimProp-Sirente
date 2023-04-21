#include "simprop.h"

using namespace simprop;

void testProductionSpectrum() {
  const double r = SI::pionMassC2 / SI::protonMassC2;
  const double eta_0 = 2. * r + r * r;
  auto numu = KelnerAharonian2008::NuMuSpectrum();
  auto antiNumu = KelnerAharonian2008::AntiNuMuSpectrum();
  auto nue = KelnerAharonian2008::NuElectronSpectrum();
  auto antiNue = KelnerAharonian2008::AntiNuElectronSpectrum();

  // KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;

  auto xAxis = utils::LogAxis<double>(1e-4, 1, 50);
  utils::OutputFile out("test_neutrino_production_spectrum.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  for (auto x : xAxis) {
    out << x << "\t";
    out << numu.Phi(1.5 * eta_0, x) / units << "\t";
    out << antiNumu.Phi(1.5 * eta_0, x) / units << "\t";
    out << nue.Phi(1.5 * eta_0, x) / units << "\t";
    out << antiNue.Phi(1.5 * eta_0, x) / units << "\t";
    out << numu.Phi(30. * eta_0, x) / units << "\t";
    out << antiNumu.Phi(30. * eta_0, x) / units << "\t";
    out << nue.Phi(30. * eta_0, x) / units << "\t";
    out << antiNue.Phi(30. * eta_0, x) / units << "\t";
    out << "\n";
  }
}

void testNeutrinoSpectrum() {
  solutions::BzNeutrinos nu;
  {
    utils::OutputFile out("SimProp_proton_spectrum.txt");
    const double units = pow2(SI::eV) / SI::m2 / SI::sr / SI::sec;
    auto E = utils::LogAxis(1e16 * SI::eV, 1e21 * SI::eV, 100);
    for (const auto& E_p : E) {
      out << std::scientific << E_p / SI::eV << "\t";
      out << pow3(E_p) * nu.getProtonFlux(E_p, 0.) / units << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("SimProp_neutrino_spectrum.txt");
    const double units = SI::GeV / SI::cm2 / SI::sr / SI::sec;
    auto E = utils::LogAxis(1e16 * SI::eV, 1e20 * SI::eV, 8 * 4);
    for (const auto& E_nu : E) {
      out << std::scientific << E_nu / SI::eV << "\t";
      out << E_nu * E_nu * nu.computeNeutrinoFlux(E_nu, 5.) / units << "\t";
      out << "\n";
    }
  }
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