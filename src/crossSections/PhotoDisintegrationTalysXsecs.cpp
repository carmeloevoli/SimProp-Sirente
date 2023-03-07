// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/crossSections/PhotoDisintegrationTalysXsecs.h"

#include "simprop/core/units.h"
#include "simprop/utils/io.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace xsecs {

TalysChannel::TalysChannel() {
  LOGD << "calling " << __func__ << " constructor";
  buildEnergyAxis();
}

void TalysChannel::buildEnergyAxis() {
  const auto minEnergy = 1. * SI::MeV;
  const auto maxEnergy = 1e2 * SI::MeV;
  const size_t sizeEnergy = 100;
  m_energyAxis = utils::LogAxis<double>(minEnergy, maxEnergy, sizeEnergy);
}

void TalysChannel::loadXsecMaps(const std::string filename) {
  auto v = utils::loadFileByRow(filename, ",");
  for (auto line : v) {
    auto pid = getPidNucleus(line[1], line[0]);
    assert(pidIsNucleus(pid));
    std::vector<double> sigma;
    std::copy(line.begin() + 2, line.end(), back_inserter(sigma));
    std::for_each(sigma.begin(), sigma.end(), [](double &x) { x *= SI::mbarn; });
    assert(sigma.size() == m_energyAxis.size());
    assert(m_xmap.insert({pid, sigma}).second);
  }
}

double TalysChannel::get(PID pid, double eps) const {
  double value = 0;
  if (!utils::isInside(eps, m_energyAxis)) return value;
  auto it = m_xmap.find(pid);
  if (it != m_xmap.end()) {
    value = utils::interpolate(eps, m_energyAxis, it->second);
  }
  return value;
}

PhotoDisintegrationTalysXsec::PhotoDisintegrationTalysXsec() {
  LOGD << "calling " << __func__ << " constructor";
  m_xsec_single.loadXsecMaps(m_singleNucleonFilename);
  m_xsec_alpha.loadXsecMaps(m_alphaFilename);
}

double PhotoDisintegrationTalysXsec::getPhotonEnergyThreshold() const { return 0.; }

double PhotoDisintegrationTalysXsec::getAtEpsPrime(PID pid, double eps) const {
  double value = m_xsec_single.get(pid, eps) + m_xsec_alpha.get(pid, eps);
  return std::max(value, 0.);
}

}  // namespace xsecs
}  // namespace simprop
