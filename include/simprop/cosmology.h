#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

#include "simprop/units.h"

namespace simprop {
namespace cosmo {

class Cosmology {
 protected:
  double m_tCmb;
  double m_h;
  double m_H0;
  double m_OmegaB;
  double m_OmegaC;
  double m_OmegaM;
  double m_OmegaL;

 public:
  Cosmology(double cmbTemperature, double littleh, double OmegaBaryon_h2, double OmegaDarkMatter_h2,
            double OmegaLambda);
  virtual ~Cosmology() = default;

  double H(double z) const;
  double hubbleTime(double z) const;
  double dtdz(double z) const;

 public:
  const double& tCmb = m_tCmb;
  const double& h = m_h;
  const double& H0 = m_H0;
  const double& OmegaB = m_OmegaB;
  const double& OmegaC = m_OmegaC;
  const double& OmegaM = m_OmegaM;
  const double& OmegaL = m_OmegaL;
};

class Planck2018 final : public Cosmology {
 public:
  Planck2018() : Cosmology(2.7255 * SI::K, 0.674, 0.02237, 0.1200, 0.685) {}
};

}  // namespace cosmo
}  // namespace simprop

#endif