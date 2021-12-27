#include "simprop/crossSections/PhotoPionProduction.h"

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

#define MAXLOGE 18.0

namespace simprop {
namespace xsecs {

double PhotoPionProduction::get(PID pid, double photonEnergy) const {
  constexpr double m_photonEnergyThreshold =
      SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  double value = 0;
  if (photonEnergy > m_photonEnergyThreshold) {
    auto loge = std::log10(photonEnergy / SI::eV);
    value = (loge < MAXLOGE) ? m_sigmas.spline(loge) : m_sigmas.spline(MAXLOGE);
  }
  return value * SI::mbarn;
}

}  // namespace xsecs
}  // namespace simprop