#include "simprop.h"

using namespace simprop;

double rate_cmb(double Gamma, double z = 0) {
  const auto cmb = photonfield::CMB();
  const auto sigma = interactions::PhotoPionProduction();
  auto value = SI::cLight / 2. / utils::pow<2>(Gamma);

  const double threshold = sigma.getPhotonEnergyThresholdCoM();
  const double epsPrimeMin = std::max(threshold, 2. * Gamma * cmb.getMinPhotonEnergy());
  const double epsPrimeMax = 2. * Gamma * cmb.getMaxPhotonEnergy();

  value *= gsl::simpsonIntegration<double>(
      [Gamma, &cmb, &sigma](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * sigma.getAtPhotonEnergyCoM(epsPrime) *
               cmb.I_gamma(epsPrime / 2. / Gamma);
      },
      std::log(epsPrimeMin), std::log(epsPrimeMax), 100);

  return value;
}

double rate_ebl(double Gamma, double z = 0) {
  const auto ebl = photonfield::Dominguez2011PhotonField();
  const auto sigma = interactions::PhotoPionProduction();
  auto value = SI::cLight / 2. / utils::pow<2>(Gamma);

  const double threshold = sigma.getPhotonEnergyThresholdCoM();
  const double epsPrimeMin = std::max(threshold, 2. * Gamma * ebl.getMinPhotonEnergy());
  const double epsPrimeMax = 2. * Gamma * ebl.getMaxPhotonEnergy();

  value *= gsl::simpsonIntegration<double>(
      [Gamma, z, &ebl, &sigma](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * sigma.getAtPhotonEnergyCoM(epsPrime) *
               ebl.I_gamma(epsPrime / 2. / Gamma, z);
      },
      std::log(epsPrimeMin), std::log(epsPrimeMax), 100);

  return value;
}

int main() {
  try {
    log::startup_information();
    const auto Gammas = utils::LogAxis(1e8, 1e16, 8 * 32);
    utils::OutputFile out("test_rates.txt");
    const auto units = 1. / SI::Mpc;
    out() << std::scientific;
    for (auto Gamma : Gammas) {
      out() << std::log10(Gamma) << "\t";
      out() << rate_cmb(Gamma) / SI::cLight / units << "\t";
      out() << rate_ebl(Gamma) / SI::cLight / units << "\t";
      out() << rate_ebl(Gamma, 3.0) / SI::cLight / units << "\t";
      out() << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
