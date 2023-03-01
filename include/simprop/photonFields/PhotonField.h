// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_PHOTONFIELDS_PHOTONFIELD_H_
#define SIMPROP_PHOTONFIELDS_PHOTONFIELD_H_

#include <memory>
#include <vector>

namespace simprop {
namespace photonfields {

class PhotonField {
 public:
  PhotonField() {}
  virtual ~PhotonField() = default;

  virtual double density(double epsRestFrame, double z = 0.) const = 0;
  virtual double I_gamma(double epsRestFrame, double z = 0.) const = 0;

  double computeIntegratedDensity(double z = 0.) const;

  virtual double getMinPhotonEnergy() const = 0;
  virtual double getMaxPhotonEnergy() const = 0;
};

using PhotonFields = std::vector<std::shared_ptr<photonfields::PhotonField>>;

}  // namespace photonfields
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H_
