#include "simprop/interactions/pppEbl.h"

#include <cmath>

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double pppEbl::rate(PID pid, double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  Gamma *= 1. + z;

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

  return pow3(1. + z) * std::max(value, 0.);
}

}  // namespace interactions
}  // namespace simprop