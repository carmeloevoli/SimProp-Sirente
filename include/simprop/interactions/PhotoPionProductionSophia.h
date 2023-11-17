#ifndef SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTIONSOPHIA_H
#define SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTIONSOPHIA_H

#include "simprop/interactions/PhotoPionProduction.h"

namespace simprop {
namespace interactions {

class PhotoPionProductionSophia final : public PhotoPionProduction {
 public:
  PhotoPionProductionSophia(const std::shared_ptr<photonfields::PhotonField>& phField);

  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;
};

}  // namespace interactions
}  // namespace simprop

#endif  // SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTION_H
