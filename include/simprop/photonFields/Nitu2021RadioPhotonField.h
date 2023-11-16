// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_PHOTONFIELDS_NITU2021RADIOPHOTONFIELD_H_
#define SIMPROP_PHOTONFIELDS_NITU2021RADIOPHOTONFIELD_H_

#include <utility>

#include "simprop/core/units.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace photonfields {

using PhotonEnergyRange = std::pair<double, double>;

class Nitu2021RadioPhotonField final : public PhotonField {
 public:
  Nitu2021RadioPhotonField();

  double density(double epsRestFrame, double z = 0.) const override;
  double I_gamma(double epsRestFrame, double z = 0.) const override;

  double getMinPhotonEnergy() const override { return m_epsRange.first; }
  double getMaxPhotonEnergy() const override { return m_epsRange.second; }

 protected:
  PhotonEnergyRange m_epsRange{0., 0.};
};

}  // namespace photonfields
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_CMBPHOTONFIELD_H_
