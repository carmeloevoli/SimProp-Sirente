#include "simprop.h"
#include "simprop/utils/misc.h"

using namespace simprop;

int main() {
  try {
    log::startup_information();
    const auto losses = losses::BGG2002ContinuousLosses();
    const auto energyAxis = utils::LogAxis(1e16 * SI::eV, 1e24 * SI::eV, 100);
    utils::OutputFile out("test_losses.txt");
    out() << std::scientific;
    for (auto E : energyAxis) {
      out() << E / SI::eV << "\t";
      out() << losses.dlnGamma_dt_0(E, proton) << "\t";
      out() << losses.dlnGamma_dz(0., E, proton) / cosmo::dtdz(0.) / (1. / SI::year) << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}