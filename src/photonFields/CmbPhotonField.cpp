#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace photonfields {

double CMB::density(double ePhoton, double z) const {
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  if (ePhoton < m_epsRange.first || ePhoton > m_epsRange.second) return 0;
  const auto kT = SI::kBoltzmann * m_temperature;
  const auto density = factor * pow2(ePhoton) / std::expm1(ePhoton / kT);
  return std::max(density, 0.);
}

double CMB::I_gamma(double ePhoton, double z) const {
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  if (ePhoton < m_epsRange.first || ePhoton > m_epsRange.second) return 0;
  const auto kT = SI::kBoltzmann * m_temperature;
  const auto I = kT * factor * (-std::log(1. - std::exp(-ePhoton / kT)));
  return std::fabs(I);
}

}  // namespace photonfields
}  // namespace simprop