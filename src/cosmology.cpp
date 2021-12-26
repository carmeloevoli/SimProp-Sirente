#include "simprop/cosmology.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace cosmo {

Cosmology::Cosmology(double tCmb, double h, double OmegaB_h2, double OmegaC_h2, double OmegaL) {
  m_tCmb = tCmb;
  m_h = h;
  m_H0 = 100. * h * SI::km / SI::sec / SI::Mpc;
  m_OmegaB = OmegaB_h2 / h / h;
  m_OmegaC = OmegaC_h2 / h / h;
  m_OmegaM = m_OmegaB + m_OmegaC;
  m_OmegaL = OmegaL;
  LOGD << "OmegaM + OmegaL = " << m_OmegaM + m_OmegaL;
  LOGD << "1/H_0 = " << 1. / m_H0 / SI::Gyr << " Gyr";
}

double Cosmology::H(double z) const {
  const auto x = 1. + z;
  return m_H0 * std::sqrt(m_OmegaM * pow3(x) + m_OmegaL);
}

double Cosmology::hubbleTime(double z) const { return 1.0 / H(z); }

double Cosmology::dtdz(double z) const { return 1. / H(z) / (1. + z); }

}  // namespace cosmo
}  // namespace simprop