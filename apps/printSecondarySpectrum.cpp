#include "simprop.h"

using namespace simprop;

const double r = SI::pionMassC2 / SI::protonMassC2;
const double eta_0 = 2. * r + r * r;

void plot_spectrum() {
  auto numu = KelnerAharonian2008::NuMuSpectrum();
  auto antiNumu = KelnerAharonian2008::AntiNuMuSpectrum();
  auto nue = KelnerAharonian2008::NuElectronSpectrum();
  auto antiNue = KelnerAharonian2008::AntiNuElectronSpectrum();
  auto xAxis = utils::LogAxis<double>(1e-4, 1, 50);
  utils::OutputFile out("test_neutrino_spectrum.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  for (auto x : xAxis) {
    out << x << "\t";
    out << numu.Phi(1.5 * eta_0, x) / units << "\t";
    out << antiNumu.Phi(1.5 * eta_0, x) / units << "\t";
    out << nue.Phi(1.5 * eta_0, x) / units << "\t";
    out << antiNue.Phi(1.5 * eta_0, x) / units << "\t";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    plot_spectrum();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}