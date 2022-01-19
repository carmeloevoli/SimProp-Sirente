#include "simprop/interactions/pppCmb.h"

#include <cmath>

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double pppCmb::rate(PID pid, double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  auto threshold = m_sigma->getThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_cmb->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_cmb->getMaxPhotonEnergy());

  value *= utils::simpsonIntegration<double>(
      [&](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * m_sigma->get(proton, epsPrime) *
               m_cmb->I_gamma(epsPrime / 2. / Gamma);
      },
      lnEpsPrimeMin, lnEpsPrimeMax, 200);

  return value;
}

}  // namespace interactions
}  // namespace simprop