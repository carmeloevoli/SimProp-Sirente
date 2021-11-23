#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace photonfield {

double CMB::getPhotonDensity(double ePhoton, double z) const {
  if (ePhoton < m_ePhotonMin || ePhoton > m_ePhotonMax) return 0;
  const auto kT = SI::kBoltzmann * m_blackbodyTemperature * (1. + z);
  const auto density = m_factor * utils::pow<2>(ePhoton) / std::expm1(ePhoton / kT);
  return density;
}

double CMB::I_gamma(double ePhoton, double z) const {
  if (ePhoton < m_ePhotonMin || ePhoton > m_ePhotonMax) return 0;
  const auto kT = SI::kBoltzmann * m_blackbodyTemperature * (1. + z);
  const auto I = kT * m_factor * (-std::log(1. - std::exp(-ePhoton / kT)));
  return std::fabs(I);
}

}  // namespace photonfield
}  // namespace simprop