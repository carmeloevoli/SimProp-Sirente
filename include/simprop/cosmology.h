#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

namespace simprop {
namespace cosmo {

#include "simprop/units.h"

class Cosmology {
 protected:
  double m_h = 0.7;
  double m_H0 = 100. * m_h * SI::km / SI::sec / SI::Mpc;
  double m_OmegaB = 0.;
  double m_OmegaC = 0.3;
  double m_OmegaM = m_OmegaC;
  double m_OmegaL = 1. - m_OmegaM;

 public:
  Cosmology(double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2, double OmegaLambda);
  virtual ~Cosmology() = default;

  double H(double z) const;
  double hubbleTime(double z) const;
  double dtdz(double z) const;

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