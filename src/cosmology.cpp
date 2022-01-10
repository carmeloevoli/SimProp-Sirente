#include "simprop/cosmology.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace cosmo {

Cosmology::Cosmology(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2,
                     double OmegaLambda) {
  m_h = littleh;
  m_H0 = 100. * m_h * SI::km / SI::sec / SI::Mpc;
  m_OmegaB = OmegaBaryon_h2 / pow2(m_h);
  m_OmegaC = OmegaDarkMatter_h2 / pow2(m_h);
  m_OmegaM = m_OmegaB + m_OmegaC;
  m_OmegaL = OmegaLambda;
  LOGD << "OmegaM + OmegaL = " << m_OmegaM + m_OmegaL;
  LOGD << "H_0^-1 = " << 1. / m_H0 / SI::Gyr << " Gyr";
}

double Cosmology::H(double z) const {
  const auto x = 1. + z;
  return m_H0 * std::sqrt(m_OmegaM * pow3(x) + m_OmegaL);
}

double Cosmology::hubbleTime(double z) const { return 1.0 / H(z); }

double Cosmology::dtdz(double z) const { return 1. / H(z) / (1. + z); }

}  // namespace cosmo
}  // namespace simprop