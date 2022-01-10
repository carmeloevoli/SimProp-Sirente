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
    m_ePhotonMax = 0.1 * SI::eV;
  }
  CMB() : CMB(2.725 * SI::K) {}

  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
  double getMinPhotonEnergy(double z = 0) const override { return m_ePhotonMin; }
  double getMaxPhotonEnergy(double z = 0) const override { return m_ePhotonMax; }

 protected:
  double m_temperature;
};

}  // namespace photonfields
}  // namespace simprop

#endif