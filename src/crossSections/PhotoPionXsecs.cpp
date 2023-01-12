#include "simprop/crossSections/PhotoPionXsecs.h"

#include "simprop/core/units.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

#define MAXLOGS 9.0

namespace simprop {
namespace xsecs {

PhotoPionXsec::PhotoPionXsec() {
  LOGD << "calling " << __func__ << " constructor";
  {
    auto filename = "data/xsecs_photopion_proton_sophia.txt";
    if (!utils::fileExists(filename)) throw std::runtime_error("data file not found");
    m_proton_sigma.loadTable(filename, 1);
    m_proton_phi.loadTable(filename, 2);
  }
  {
    auto filename = "data/xsecs_photopion_neutron_sophia.txt";
    if (!utils::fileExists(filename)) throw std::runtime_error("data file not found");
    m_neutron_sigma.loadTable(filename, 1);
    m_neutron_phi.loadTable(filename, 2);
  }
}

double PhotoPionXsec::getPhotonEnergyThreshold() const {
  constexpr double photonEnergyThreshold =
      SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  return photonEnergyThreshold;
}

double PhotoPionXsec::getProtonXsec(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s > sThreshold && m_proton_sigma.xIsInside(s / SI::GeV2)) {
    value = m_proton_sigma.get(s / SI::GeV2);
  }
  return std::max(value, 0.) * SI::mbarn;
}

double PhotoPionXsec::getNeutronXsec(double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::neutronMassC2 + SI::pionMassC2);
  if (s > sThreshold && m_neutron_sigma.xIsInside(s / SI::GeV2)) {
    value = m_neutron_sigma.get(s / SI::GeV2);
  }
  return std::max(value, 0.) * SI::mbarn;
}

double PhotoPionXsec::getAtS(PID pid, double s) const {
  auto Z = getPidNucleusCharge(pid);
  auto A = getPidNucleusMassNumber(pid);
  auto value = (double)Z * getProtonXsec(s);
  if (A > Z) value += (double)(A - Z) * getNeutronXsec(s);
  return value;
}

double PhotoPionXsec::getPhiAtS(PID pid, double s) const {
  double value = 0;
  constexpr auto sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (s < sThreshold) return value;
  if (pid == proton && m_proton_phi.xIsInside(s / SI::GeV2)) {
    value = m_proton_phi.get(s / SI::GeV2);
  }
  if (pid == neutron && m_neutron_phi.xIsInside(s / SI::GeV2)) {
    value = m_neutron_phi.get(s / SI::GeV2);
  }
  if (!pidIsNucleon(pid)) throw std::runtime_error("phi not implemented for nuclei in the SPM");
  return std::max(value, 0.) * pow2(SI::GeV2) * SI::mbarn;
}

}  // namespace xsecs
}  // namespace simprop