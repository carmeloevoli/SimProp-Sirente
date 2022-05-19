#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

#include <vector>

#include "simprop/units.h"

namespace simprop {
namespace cosmo {

class Cosmology {
 protected:
  double m_h = 0.7;
  double m_H0 = 70. * m_h * SI::km / SI::sec / SI::Mpc;
  double m_OmegaB = 0.;
  double m_OmegaC = 0.3;
  double m_OmegaM = m_OmegaC;
  double m_OmegaL = 1. - m_OmegaM;
  double m_Dh = SI::cLight / m_H0;
  double m_age = -1;

  const size_t m_size = 1000;
  const double m_zmin = 0.0001;
  const double m_zmax = 100;

  std::vector<double> m_z;   // redshift
  std::vector<double> m_Dc;  // comoving distance
  std::vector<double> m_Dl;  // luminosity distance
  std::vector<double> m_Dt;  // light travel distance

 public:
  Cosmology();
  Cosmology(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2, double OmegaLambda);
  virtual ~Cosmology() = default;

  inline double E(double z) const { return std::sqrt(m_OmegaL + m_OmegaM * pow3(1. + z)); }
  inline double hubbleRate(double z) const { return m_H0 * E(z); }
  inline double dtdz(double z) const { return 1. / hubbleRate(z) / (1. + z); }

  double comovingDistance2Redshift(double d) const;
  double redshift2ComovingDistance(double z) const;
  double luminosityDistance2Redshift(double d) const;
  double redshift2LuminosityDistance(double z) const;
  double lightTravelDistance2Redshift(double d) const;
  double redshift2LightTravelDistance(double z) const;
  double comoving2LightTravelDistance(double d) const;
  double lightTravel2ComovingDistance(double d) const;
  double redshift2UniverseAge(double z) const;

 protected:
  void setParameters(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2,
                     double OmegaLambda);
  void init();
  double computeLightTravelDistance(double z) const;
  double computeComovingDistance(double z) const;

 public:
  const double& h = m_h;
  const double& H0 = m_H0;
  const double& OmegaB = m_OmegaB;
  const double& OmegaC = m_OmegaC;
  const double& OmegaM = m_OmegaM;
  const double& OmegaL = m_OmegaL;
};

class Planck2018 final : public Cosmology {
 public:
  Planck2018() : Cosmology(0.674, 0.02237, 0.1200, 0.685) {}
};

}  // namespace cosmo
}  // namespace simprop

#endif