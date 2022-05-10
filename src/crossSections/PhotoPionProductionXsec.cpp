#include "simprop/crossSections/PhotoPionProductionXsec.h"

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

#define MAXLOGS 9.0

namespace simprop {
namespace xsecs {

PhotoPionProductionXsec::PhotoPionProductionXsec() {
  LOGD << "calling " << __func__ << " constructor";
  m_sigma.loadTable(m_filename, 1);
  m_phi.loadTable(m_filename, 2);
}

double PhotoPionProductionXsec::getPhotonEnergyThreshold() const {
  constexpr double photonEnergyThreshold =
      SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  return photonEnergyThreshold;
}

double PhotoPionProductionXsec::getAtS(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold && m_sigma.xIsInside(s / SI::GeV2)) {
    value = m_sigma.get(s / SI::GeV2);
  }
  return std::max(value, 0.) * SI::mbarn;
}

double PhotoPionProductionXsec::getAtEpsPrime(double epsPrime) const {
  if (epsPrime > getPhotonEnergyThreshold()) {
    auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
    return getAtS(s);
  } else {
    return 0;
  }
}

double PhotoPionProductionXsec::getPhiAtS(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold && m_phi.xIsInside(s / SI::GeV2)) {
    value = m_phi.get(s / SI::GeV2);
  }
  return std::max(value, 0.) * pow2(SI::GeV2) * SI::mbarn;
}

}  // namespace xsecs
}  // namespace simprop