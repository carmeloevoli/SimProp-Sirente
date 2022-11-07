#ifndef SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_PHOTONFIELD_H

#include <memory>
#include <vector>

namespace simprop {
namespace photonfields {

class PhotonField {
 public:
  PhotonField() {}
  virtual ~PhotonField() = default;

  virtual double density(double ePhoton, double z = 0.) const = 0;
  virtual double I_gamma(double ePhoton, double z = 0.) const = 0;

  virtual double getMinPhotonEnergy() const = 0;
  virtual double getMaxPhotonEnergy() const = 0;
};

using PhotonFields = std::vector<std::shared_ptr<photonfields::PhotonField>>;

}  // namespace photonfields
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
