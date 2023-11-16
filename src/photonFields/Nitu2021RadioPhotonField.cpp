// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/photonFields/Nitu2021RadioPhotonField.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace photonfields {

Nitu2021RadioPhotonField::Nitu2021RadioPhotonField() {
  m_epsRange.first = 1e3 * SI::Hz * SI::hPlanck;  // Hertz
  m_epsRange.second = 1e12 * SI::Hz * SI::hPlanck;
  LOGD << "calling " << __func__ << " constructor";
}

double Nitu2021RadioPhotonField::density(double epsRestFrame, double z) const {
  if (epsRestFrame < m_epsRange.first || epsRestFrame > m_epsRange.second) return 0;
  // Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
  std::vector<double> p;
  p.push_back(-1.9847e1);
  p.push_back(-2.9857e-1);
  p.push_back(-2.6984e-1);
  p.push_back(9.5394e-2);
  p.push_back(-4.9059e-2);
  p.push_back(4.4297e-3);
  p.push_back(7.6038e-3);
  p.push_back(-1.9690e-3);
  p.push_back(-2.2573e-4);
  p.push_back(1.1762e-4);
  p.push_back(-9.9443e-6);

  auto nu = epsRestFrame / SI::hPlanck;
  auto I = 0.;
  for (size_t k = 0; k < p.size(); ++k) {
    I += p[k] * std::pow(std::log10(nu / SI::MHz), k);
  }
  I = std::pow(10., I) * SI::watt / SI::Hz / pow2(SI::m2) / SI::sr;
  I = 4 * M_PI / (SI::hPlanck * SI::cLight) * (I / epsRestFrame);

  return I * std::pow(1. + z, 3.);
}

double Nitu2021RadioPhotonField::I_gamma(double epsRestFrame, double z) const {
  if (epsRestFrame < m_epsRange.first || epsRestFrame > m_epsRange.second) return 0;
  return 0;  // TODO to be implemented
}

}  // namespace photonfields
}  // namespace simprop
