#include "simprop/crossSections/PhotoPionProductionXsec.h"

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

#define MAXLOGS 9.0

namespace simprop {
namespace xsecs {

double PhotoPionProductionXsec::getPhotonEnergyThreshold() const {
  constexpr double photonEnergyThreshold =
      SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  return photonEnergyThreshold;
}

double PhotoPionProductionXsec::getAtS(PID pid, double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold) {
    auto logs = std::log10(s / pow2(SI::GeV));
    value = (logs < MAXLOGS) ? m_sigmas.get(logs) : m_sigmas.get(MAXLOGS);
  }
  return std::max(value, 0.) * SI::mbarn;
}

double PhotoPionProductionXsec::getAtEpsPrime(PID pid, double epsPrime) const {
  if (epsPrime > getPhotonEnergyThreshold()) {
    auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
    return getAtS(pid, s);
  } else {
    return 0;
  }
}

}  // namespace xsecs
}  // namespace simprop