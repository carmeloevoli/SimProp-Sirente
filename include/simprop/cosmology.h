#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

#include "simprop/units.h"

namespace simprop {
namespace cosmo {
/**
 * \addtogroup Cosmology
 * @{
 */

/**
 @class Cosmology cosmology.h include/simprop/cosmology.h
 @brief It contains all the functions related to the Universe evolution.
 */
class Cosmology {
 protected:
  double m_tCmb = 2.7 * SI::K;                         /**< CMB Temperature at \f$z=0\f$ */
  double m_h = 0.7;                                    /**< reduced Hubble constant */
  double m_H0 = 100. * h * SI::km / SI::sec * SI::Mpc; /**< Hubble constant */
  double m_OmegaB = 0.;                                /**< Baryon density parameter */
  double m_OmegaC = 0.3;                               /**< Dark matter density parameter */
  double m_OmegaM = 0.7;                               /**< Matter density parameter */
  double m_OmegaL = 0.3;                               /**< Dark energy density parameter */

 public:
  /**
   * @brief Construct the Cosmology object using 5 cosmological parameter
   *
   * @param tCmb CMB temperature at \f$z = 0\f$
   * @param h reduced Hubble constant
   * @param OmegaB_h2 \f$\Omega_b / h^2\f$
   * @param OmegaC_h2 \f$\Omega_c / h^2\f$
   * @param OmegaL \f$\Omega_\Lambda\f$
   */
  Cosmology(double tCmb, double h, double OmegaB_h2, double OmegaC_h2, double OmegaL);

  /**
   * @brief Destroy the Cosmology object
   */
  virtual ~Cosmology() = default;

  /**
   * @brief compute the Hubble "constant" at z
   *
   * @param z redshift
   * @return double \f$H(z)\f$
   */
  double H(double z) const;

  /**
   * @brief returns hubble time \f$t_h = 1/H\f$
   *
   * @param z redsfhit
   */
  double hubbleTime(double z) const;

  /**
   * @brief returns the value of \f$dt/dz\f$ at the redshift parameter z
   *
   * @param z redsfhit
   */
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
  /**
   * @brief Construct the Planck2018 Cosmology object using
   * cosmological parameter from REF.
   */
  Planck2018() : Cosmology(2.7255 * SI::K, 0.674, 0.02237, 0.1200, 0.685) {}
};

}  // namespace cosmo
}  // namespace simprop

#endif