#include "simprop/photonFields/CmbPhotonField.h"

#include "simprop/utils/io.h"

namespace simprop {
namespace photonfields {

CMB::CMB(double T)
    : m_temperature(T),
      m_epsRange({1e-5 * (T / 2.725 / SI::K) * SI::eV, 0.1 * (T / 2.725 / SI::K) * SI::eV}) {
  LOGD << "calling " << __func__ << " constructor";
}

double CMB::density(double ePhoton, double z) const {
  if (ePhoton < m_epsRange.first || ePhoton > m_epsRange.second) return 0;
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  const auto kT = SI::kBoltzmann * m_temperature;
  const auto density = factor * pow2(ePhoton) / std::expm1(ePhoton / kT);
  return std::max(density, 0.);
}

double CMB::I_gamma(double ePhoton, double z) const {
  if (ePhoton < m_epsRange.first || ePhoton > m_epsRange.second) return 0;
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  const auto kT = SI::kBoltzmann * m_temperature;
  const auto I = kT * factor * (-std::log(1. - std::exp(-ePhoton / kT)));
  return std::fabs(I);
}

}  // namespace photonfields
}  // namespace simprop