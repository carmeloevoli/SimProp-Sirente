#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace photonfields {

double CMB::density(double ePhoton, double z) const {
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  if (ePhoton < m_ePhotonMin || ePhoton > m_ePhotonMax) return 0;
  const auto kT = SI::kBoltzmann * m_temperature * (1. + z);
  const auto density = factor * pow2(ePhoton) / std::expm1(ePhoton / kT);
  return density * pow3(1. + z);
}

double CMB::I_gamma(double ePhoton, double z) const {
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  if (ePhoton < m_ePhotonMin || ePhoton > m_ePhotonMax) return 0;
  const auto kT = SI::kBoltzmann * m_temperature * (1. + z);
  const auto I = kT * factor * (-std::log(1. - std::exp(-ePhoton / kT)));
  return std::fabs(I) * pow3(1. + z);
}

}  // namespace photonfields
}  // namespace simprop