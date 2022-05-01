#include "simprop/cosmology.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace cosmo {

Cosmology::Cosmology() { init(); }

Cosmology::Cosmology(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2,
                     double OmegaLambda) {
  setParameters(littleh, OmegaBaryon_h2, OmegaDarkMatter_h2, OmegaLambda);
  init();
  LOGD << "calling " << __func__ << " constructor";
  LOGD << "OmegaM + OmegaL = " << m_OmegaM + m_OmegaL;
  LOGD << "H_0^-1 = " << 1. / m_H0 / SI::Gyr << " Gyr";
  LOGD << "d_H = " << m_Dh / SI::Mpc << " Mpc";
}

void Cosmology::setParameters(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2,
                              double OmegaLambda) {
  m_h = littleh;
  m_H0 = 100. * m_h * SI::km / SI::sec / SI::Mpc;
  m_OmegaB = OmegaBaryon_h2 / pow2(m_h);
  m_OmegaC = OmegaDarkMatter_h2 / pow2(m_h);
  m_OmegaM = m_OmegaB + m_OmegaC;
  m_OmegaL = OmegaLambda;
  m_Dh = SI::cLight / m_H0;
}

void Cosmology::init() {
  m_z.resize(m_size);
  m_Dc.resize(m_size);
  m_Dl.resize(m_size);
  m_Dt.resize(m_size);

  m_z[0] = 0;
  m_Dc[0] = 0;
  m_Dl[0] = 0;
  m_Dt[0] = 0;

  const double dlz = std::log10(m_zmax) - std::log10(m_zmin);
  for (size_t i = 1; i < m_size; i++) {
    const double z = m_zmin * std::pow(10., i * dlz / (m_size - 1));
    m_z[i] = z;
    m_Dc[i] = computeComovingDistance(z);
    m_Dl[i] = (1. + z) * m_Dc[i];
    m_Dt[i] = computeLightTravelDistance(z);
  }
}

double Cosmology::computeLightTravelDistance(double z) const {
  auto integrand = [&](double x) { return 1. / E(x) / (1. + x); };
  return m_Dh * utils::QAGIntegration<double>(integrand, 0., z, 1000, 1e-8);
}

double Cosmology::computeComovingDistance(double z) const {
  auto integrand = [&](double x) { return 1. / E(x); };
  return m_Dh * utils::QAGIntegration<double>(integrand, 0., z, 1000, 1e-8);
}

double Cosmology::comovingDistance2Redshift(double d) const {
  if (d < 0) throw std::invalid_argument("Cosmology: d < 0");
  if (d > m_Dc.back()) throw std::invalid_argument("Cosmology: d > dmax");
  return utils::interpolate(d, m_Dc, m_z);
}

double Cosmology::redshift2ComovingDistance(double z) const {
  if (z < 0) throw std::invalid_argument("Cosmology: z < 0");
  if (z > m_zmax) throw std::invalid_argument("Cosmology: z > zmax");
  return utils::interpolate(z, m_z, m_Dc);
}

double Cosmology::luminosityDistance2Redshift(double d) const {
  if (d < 0) throw std::invalid_argument("Cosmology: d < 0");
  if (d > m_Dl.back()) throw std::invalid_argument("Cosmology: d > dmax");
  return utils::interpolate(d, m_Dl, m_z);
}

double Cosmology::redshift2LuminosityDistance(double z) const {
  if (z < 0) throw std::invalid_argument("Cosmology: z < 0");
  if (z > m_zmax) throw std::invalid_argument("Cosmology: z > zmax");
  return utils::interpolate(z, m_z, m_Dl);
}

double Cosmology::lightTravelDistance2Redshift(double d) const {
  if (d < 0) throw std::invalid_argument("Cosmology: d < 0");
  if (d > m_Dt.back()) throw std::invalid_argument("Cosmology: d > dmax");
  return utils::interpolate(d, m_Dt, m_z);
}

double Cosmology::redshift2LightTravelDistance(double z) const {
  if (z < 0) throw std::invalid_argument("Cosmology: z < 0");
  if (z > m_zmax) throw std::invalid_argument("Cosmology: z > zmax");
  return utils::interpolate(z, m_z, m_Dt);
}

double Cosmology::comoving2LightTravelDistance(double d) const {
  if (d < 0) throw std::invalid_argument("Cosmology: d < 0");
  if (d > m_Dc.back()) throw std::invalid_argument("Cosmology: d > dmax");
  return utils::interpolate(d, m_Dc, m_Dt);
}

double Cosmology::lightTravel2ComovingDistance(double d) const {
  if (d < 0) throw std::invalid_argument("Cosmology: d < 0");
  if (d > m_Dt.back()) throw std::invalid_argument("Cosmology: d > dmax");
  return utils::interpolate(d, m_Dt, m_Dc);
}

}  // namespace cosmo
}  // namespace simprop