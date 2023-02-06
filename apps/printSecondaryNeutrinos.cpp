#include "simprop.h"

using namespace simprop;

void plot_production_spectrum() {
  const double r = SI::pionMassC2 / SI::protonMassC2;
  const double eta_0 = 2. * r + r * r;
  auto numu = KelnerAharonian2008::NuMuSpectrum();
  auto antiNumu = KelnerAharonian2008::AntiNuMuSpectrum();
  auto nue = KelnerAharonian2008::NuElectronSpectrum();
  auto antiNue = KelnerAharonian2008::AntiNuElectronSpectrum();
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

void test_neutrino_spectrum() {
  solutions::Beniamino b(true);
  b.setSlope(2.7);
  b.setSourceEvolution(0.);

  double I = 0;
  for (int i = 0; i < 10000; ++i) {
    I = b.I_dEpsilon(1e17 * SI::eV, 1e19 * SI::eV, 0.1);
  }
  std::cout << I << "\n";
}

void plot_neutrino_spectrum() {
  // auto eAxis = utils::LogAxis<double>(1e-4, 1, 50);
  // utils::OutputFile out("test_neutrino_observed_spectrum.txt");
  // out << "# x - spectrum\n";
  // out << std::scientific;
  solutions::Beniamino b(true);
  b.setSlope(2.7);
  b.setSourceEvolution(0.);
  auto E = 1e8 * SI::GeV;
  const auto units = SI::GeV / SI::cm2 / SI::sec / SI::sr;
  LOGD << E / SI::GeV << " " << E * E * b.computeNeutrinoFlux(1e8 * SI::GeV, 1.) / units;
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    // plot_production_spectrum();
    test_neutrino_spectrum();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}