#ifndef SIMPROP_PHOTONFIELDS_CMB_H
#define SIMPROP_PHOTONFIELDS_CMB_H

#include <cmath>
#include <string>

#include "simprop/photonFields/AbstractPhotonField.h"
#include "simprop/units.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace photonfield {

class CMB : public AbstractPhotonField {
 public:
  CMB() : AbstractPhotonField("CMB", 1e-10 * SI::eV, 1e-1 * SI::eV) {}
  double getPhotonDensity(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;

 protected:
  const double m_blackbodyTemperature = 2.73 * SI::K;
  const double m_factor = 1. / utils::pow<2>(M_PI) / utils::pow<3>(SI::hbarC);
};

}  // namespace photonfield
}  // namespace simprop

#endif