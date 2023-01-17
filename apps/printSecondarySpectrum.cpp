#include "simprop.h"

using namespace simprop;

const double r = SI::pionMassC2 / SI::protonMassC2;
const double eta_0 = 2. * r + r * r;

void plot_spectrum() {
  auto spectrum = KelnerAharonian2008::AntiNuMuSecondarySpectrum();
  auto xAxis = utils::LogAxis<double>(1e-4, 1, 50);
  utils::OutputFile out("test_antinumu_spectrum.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  for (auto x : xAxis) {
    out << x << "\t";
    out << spectrum.Phi(1.5 * eta_0, x) / units << "\t";
    out << spectrum.Phi(30. * eta_0, x) / units << "\t";
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