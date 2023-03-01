// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/photonFields/CmbPhotonField.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace photonfields {

CMB::CMB(double T) : m_temperature(T) {
  m_epsRange.first = 1e-5 * SI::eV * (T / 2.725 / SI::K);
  m_epsRange.second = 0.1 * SI::eV * (T / 2.725 / SI::K);
  LOGD << "calling " << __func__ << " constructor";
}

double CMB::density(double epsRestFrame, double z) const {
  if (epsRestFrame < m_epsRange.first || epsRestFrame > m_epsRange.second) return 0;
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  const auto kT = SI::kBoltzmann * m_temperature * (1. + z);
  const auto density = factor * pow2(epsRestFrame) / std::expm1(epsRestFrame / kT);
  return std::max(density, 0.);
}

double CMB::I_gamma(double epsRestFrame, double z) const {
  if (epsRestFrame < m_epsRange.first || epsRestFrame > m_epsRange.second) return 0;
  constexpr double factor = 1. / pow2(M_PI) / pow3(SI::hbarC);
  const auto kT = SI::kBoltzmann * m_temperature * (1. + z);
  const auto I = factor * kT * (-std::log(1. - std::exp(-epsRestFrame / kT)));
  return std::fabs(I);
}

}  // namespace photonfields
}  // namespace simprop
