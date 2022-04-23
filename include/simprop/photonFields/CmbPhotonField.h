#ifndef SIMPROP_PHOTONFIELDS_CMB_H
#define SIMPROP_PHOTONFIELDS_CMB_H

#include <cmath>
#include <string>

#include "simprop/photonFields/PhotonField.h"
#include "simprop/units.h"

namespace simprop {
namespace photonfields {

using PhotonEnergyRange = std::pair<double, double>;

class CMB final : public PhotonField {
 public:
  CMB(double T);
  CMB() : CMB(2.725 * SI::K) {}

  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
  double getMinPhotonEnergy(double z = 0) const override { return m_epsRange.first; }
  double getMaxPhotonEnergy(double z = 0) const override { return m_epsRange.second; }

 protected:
  double m_temperature;
  PhotonEnergyRange m_epsRange{0., 0.};
};

}  // namespace photonfields
}  // namespace simprop

#endif