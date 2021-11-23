#include "simprop/photonFields/Dominguez2011PhotonField.h"

#include "simprop/Units.h"
#include "simprop/common.h"
#include "simprop/utils/gsl.h"

namespace simprop {
namespace photonfield {

Dominguez2011PhotonField::Dominguez2011PhotonField()
    : AbstractPhotonField("Dominguez2011", 1.240e-3 * SI::eV, 1.227 * SI::eV) {}

double Dominguez2011PhotonField::getPhotonDensity(double ePhoton, double z) const {
  const auto loge = std::log10(ePhoton / SI::eV);
  if (m_field.isWithinXRange(loge) && m_field.isWithinYRange(z)) {
    const auto value = m_field.get(z, loge);
    const auto power = std::pow(10., value) / SI::eV / SI::m3;
    return power;
  } else {
    return 0;
  }
};

double Dominguez2011PhotonField::I_gamma(double ePhoton, double z) const {
  if (ePhoton > m_ePhotonMax) return 0;
  const auto epsUpper = m_ePhotonMax;
  const auto epsLower = std::max(ePhoton, m_ePhotonMin);
  auto I = gsl::QAGIntegration<double>(
      [&](double eps) { return getPhotonDensity(eps, z) / utils::pow<2>(eps); }, epsLower, epsUpper,
      1000, 1e-3);
  return I;
}

}  // namespace photonfield
}  // namespace simprop