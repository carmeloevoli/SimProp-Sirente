#ifndef SIMPROP_PHOTONFIELDS_CMB_H
#define SIMPROP_PHOTONFIELDS_CMB_H

#include <cmath>
#include <string>

#include "simprop/photonFields/PhotonField.h"
#include "simprop/units.h"

namespace simprop {
namespace photonfields {

class CMB final : public PhotonField {
 public:
  CMB(double T) : m_temperature(T) {
    m_ePhotonMin = 1e-10 * SI::eV;
    m_ePhotonMax = 1e-1 * SI::eV;
  }
  CMB() {}

  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;

 protected:
  double m_temperature = 2.725 * SI::K;
};

}  // namespace photonfields
}  // namespace simprop

#endif