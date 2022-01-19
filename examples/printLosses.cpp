#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto cosmology = std::make_shared<cosmo::Planck2018>();
    auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
    auto losses = losses::BGG2002ContinuousLosses(cosmology);

    std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{
        std::make_shared<photonfields::CMB>(),
        std::make_shared<photonfields::Dominguez2011PhotonField>()};

    auto pair = losses::PairProductionLosses(cosmology, phFields);

    auto energyAxis = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 500);
    utils::OutputFile out("test_losses.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / SI::eV << "\t";
      auto Gamma = E / SI::protonMassC2;
      out << SI::cLight / adiabatic.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << SI::cLight / losses.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << SI::cLight / pair.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}