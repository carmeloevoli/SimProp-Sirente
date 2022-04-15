#include "simprop/crossSections/PhotoPionProductionXsec.h"

#include "simprop/units.h"
#include "simprop/utils/numeric.h"

#define MAXLOGS 9.0

namespace simprop {
namespace xsecs {

PhotoPionProductionXsec::PhotoPionProductionXsec() { loadDataFile(); }

void PhotoPionProductionXsec::loadDataFile() {
  auto v = utils::loadFileByRow(m_filename, ",");
  size_t counter = 0;
  for (size_t i = 0; i < m_sSize; ++i) {
    auto line = v.at(counter);
    if (line.size() != 3) throw std::runtime_error("error in reading table values");
    m_sEnergies.emplace_back(line[0] * SI::GeV2);
    m_sigma.emplace_back(line[1] * SI::mbarn);
    m_phi.emplace_back(line[2] * SI::mbarn * pow4(SI::GeV));
    counter++;
  }
  assert(m_sEnergies.size() == m_sSize);
}

double PhotoPionProductionXsec::getPhotonEnergyThreshold() const {
  constexpr double photonEnergyThreshold =
      SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  return photonEnergyThreshold;
}

double PhotoPionProductionXsec::getAtS(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold && utils::isInside(s, m_sEnergies)) {
    value = utils::interpolate(s, m_sEnergies, m_sigma);
  }
  return std::max(value, 0.);
}

double PhotoPionProductionXsec::getPhiAtS(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold && utils::isInside(s, m_sEnergies)) {
    value = utils::interpolate(s, m_sEnergies, m_phi);
  }
  return std::max(value, 0.);
}

double PhotoPionProductionXsec::getAtEpsPrime(double epsPrime) const {
  if (epsPrime > getPhotonEnergyThreshold()) {
    auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
    return getAtS(s);
  } else {
    return 0;
  }
}

}  // namespace xsecs
}  // namespace simprop