// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/crossSections/PhotoDisintegrationTalysXsecs.h"

#include "simprop/core/units.h"
#include "simprop/utils/logging.h"

namespace simprop {
namespace xsecs {

PhotoDisintegrationTalysXsec::PhotoDisintegrationTalysXsec() {
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoDisintegrationTalysXsec::getPhotonEnergyThreshold() const { return 0.; }

double PhotoDisintegrationTalysXsec::getAtS(PID pid, double s) const {
  double value = 0;
  return std::max(value, 0.) * SI::mbarn;
}

// double PhotoDisintegrationTalysXsec::getAtEpsPrime(double epsPrime) const {
//   if (epsPrime > getPhotonEnergyThreshold()) {
//     auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
//     return getAtS(s);
//   } else {
//     return 0;
//   }
// }

double PhotoDisintegrationTalysXsec::getPhiAtS(PID pid, double s) const {
  double value = 0;
  return std::max(value, 0.) * pow2(SI::GeV2) * SI::mbarn;
}

}  // namespace xsecs
}  // namespace simprop
