#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto cosmology = std::make_shared<cosmo::Planck2018>();
    auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
    auto losses = losses::BGG2002ContinuousLosses();

    // std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{
    //     std::make_shared<photonfields::CMB>(),
    //     std::make_shared<photonfields::Dominguez2011PhotonField>()};

    // auto pair = losses::PairProductionLosses(cosmology, phFields);

    std::vector<std::shared_ptr<photonfields::PhotonField> > cmb{
        std::make_shared<photonfields::CMB>()};
    auto pair_cmb = losses::PairProductionLosses(cmb);

    std::vector<std::shared_ptr<photonfields::PhotonField> > ebl{
        std::make_shared<photonfields::Dominguez2011PhotonField>()};
    auto pair_irb = losses::PairProductionLosses(ebl);

    const auto gammaAxis = utils::LogAxis(1e8, 1e16, 8 * 32);

    utils::OutputFile out("test_losses.txt");
    out << std::scientific;
    for (auto Gamma : gammaAxis) {
      out << Gamma << "\t";
      out << SI::cLight / adiabatic.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << SI::cLight / losses.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << SI::cLight / pair_cmb.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << SI::cLight / pair_irb.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}