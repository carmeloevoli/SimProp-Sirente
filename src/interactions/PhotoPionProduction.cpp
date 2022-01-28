#include "simprop/interactions/PhotoPionProduction.h"

#include <cmath>

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double PhotoPionProduction::computeRateComoving(PID pid, double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  auto threshold = m_sigma->getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_ebl->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_ebl->getMaxPhotonEnergy());

  value *= utils::simpsonIntegration<double>(
      [&](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * m_sigma->getAtEpsPrime(proton, epsPrime) *
               m_ebl->I_gamma(epsPrime / 2. / Gamma, z);
      },
      lnEpsPrimeMin, lnEpsPrimeMax, 300);

  return std::max(value, 0.);
}

double PhotoPionProduction::rate(PID pid, double Gamma, double z) const {
  return pow3(1. + z) * computeRateComoving(pid, Gamma * (1. + z), z);
}

double PhotoPionProduction::sampleS(double r, double sMax) const {
  const double sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  // inline double sMax(double protonEnergy, double photonEnergy) const {
  //   return pow2(SI::protonMassC2) + 4. * protonEnergy * photonEnergy;
  // }
  // auto I_max = utils::QAGIntegration<double>(
  //     [&](double s) { return (s - pow2(SI::protonMassC2)) * m_sigma->getAtS(s); }, sThreshold,
  //     sMax, 1000, 1e-5);
  return 0;
}

}  // namespace interactions
}  // namespace simprop