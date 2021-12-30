#include "simprop.h"

using namespace simprop;

double rate_cmb(double Gamma, double z = 0) {
  const auto cmb = photonfields::CMB();
  const auto sigma = xsecs::PhotoPionProduction();
  const auto epsPrimeMin = 2. * Gamma * cmb.getMinPhotonEnergy();
  const auto epsPrimeMax = 2. * Gamma * cmb.getMaxPhotonEnergy();

  auto value = 1. / 2. / pow2(Gamma);
  auto I = utils::simpsonIntegration<double>(
      [Gamma, &cmb, &sigma](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * sigma.get(proton, epsPrime) *
               cmb.I_gamma(epsPrime / 2. / Gamma);
      },
      std::log(epsPrimeMin), std::log(epsPrimeMax), 100);
  return value * I;
}

int main() {
  try {
    utils::startup_information();
    auto cosmology = cosmo::Planck2018();
    auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
    auto losses = losses::BGG2002ContinuousLosses(cosmology);
    auto energyAxis = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 500);
    utils::OutputFile out("test_losses.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / SI::eV << "\t";
      out << SI::cLight / adiabatic.dlnE_dt(proton, E) / SI::Mpc << "\t";
      out << SI::cLight / losses.dlnE_dt(proton, E) / SI::Mpc << "\t";
      out << SI::cLight / rate_cmb(E / SI::protonMassC2) / SI::Mpc << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}