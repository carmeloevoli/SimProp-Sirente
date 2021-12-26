#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

#include "simprop/units.h"

namespace simprop {
namespace cosmo {

class Cosmology {
 protected:
  double m_tCmb = 0;
  double m_h = 0;
  double m_H0 = 0;
  double m_OmegaB = 0;
  double m_OmegaC = 0;
  double m_OmegaM = 0;
  double m_OmegaL = 0;

 public:
  Cosmology(double tCmb, double h, double OmegaB_h2, double OmegaC_h2, double OmegaL);
  virtual ~Cosmology() = default;

  /* returns the hubble "constant" at z */
  double H(double z) const;

  /* returns hubble time, t_h = 1/H */
  double hubbleTime(double z) const;

  /* returns the value of dt/dz at the redshift parameter z. */
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