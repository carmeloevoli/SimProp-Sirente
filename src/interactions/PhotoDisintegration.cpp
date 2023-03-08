#include "simprop/interactions/PhotoDisintegration.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

PhotoDisintegration::PhotoDisintegration(const std::shared_ptr<photonfields::PhotonField>& phField)
    : Interaction(phField) {
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoDisintegration::interactionLength(PID pid, double Gamma) const {
  return SI::cLight * getPidNucleusMassNumber(pid) / rate(pid, Gamma);
}

double PhotoDisintegration::rate(PID pid, double Gamma, double z) const {
  auto value = double(0);
  auto threshold = m_xs.getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_phField->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_phField->getMaxPhotonEnergy());
  if (lnEpsPrimeMax > lnEpsPrimeMin) {
    value = utils::simpsonIntegration<double>(
        [this, pid, Gamma, z](double lnEpsPrime) {
          auto epsPrime = std::exp(lnEpsPrime);
          return epsPrime * epsPrime * m_xs.getAtEpsPrime(pid, epsPrime) *
                 m_phField->I_gamma(epsPrime / 2. / Gamma, z);
        },
        lnEpsPrimeMin, lnEpsPrimeMax, 300);
    value *= SI::cLight / 2. / pow2(Gamma);
  }
  return std::max(value, 0.);
}

std::vector<Particle> PhotoDisintegration::finalState(const Particle& particle,
                                                      double zInteractionPoint,
                                                      RandomNumberGenerator& rng) const {
  return std::vector<Particle>();
}

}  // namespace interactions
}  // namespace simprop
