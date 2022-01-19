#ifndef SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_PHOTONFIELD_H

namespace simprop {
namespace photonfields {

class PhotonField {
 public:
  PhotonField() {}
  virtual ~PhotonField() = default;

  virtual double density(double ePhoton, double z = 0.) const = 0;
  virtual double I_gamma(double ePhoton, double z = 0.) const = 0;

  virtual double getMinPhotonEnergy(double z = 0) const = 0;
  virtual double getMaxPhotonEnergy(double z = 0) const = 0;
};

}  // namespace photonfields
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
