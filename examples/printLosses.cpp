#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto cosmology = cosmo::Planck2018();
    auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
    auto losses = losses::BGG2002ContinuousLosses(cosmology);
    auto pair = losses::PairProductionLosses(cosmology);
    auto energyAxis = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 500);
    utils::OutputFile out("test_losses.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / SI::eV << "\t";
      out << SI::cLight / adiabatic.dlnE_dt(proton, E) / SI::Mpc << "\t";
      out << SI::cLight / losses.dlnE_dt(proton, E) / SI::Mpc << "\t";
      out << SI::cLight / pair.dlnE_dt(proton, E) / SI::Mpc << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}