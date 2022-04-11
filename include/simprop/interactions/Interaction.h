#ifndef SIMPROP_INTERACTIONS_INTERACTION_H
#define SIMPROP_INTERACTIONS_INTERACTION_H

#include <memory>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/particle.h"
#include "simprop/photonFields/PhotonField.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

class Interaction {
 protected:
  std::shared_ptr<xsecs::CrossSection> m_sigma;
  std::shared_ptr<photonfields::PhotonField> m_phField;

 public:
  Interaction(const std::shared_ptr<xsecs::CrossSection>& sigma,
              const std::shared_ptr<photonfields::PhotonField>& phField)
      : m_sigma(sigma), m_phField(phField) {}
  virtual ~Interaction() = default;
  virtual double rate(PID pid, double Gamma, double z = 0) const = 0;
  virtual std::vector<Particle> finalState(const Particle& incomingParticle,
                                           double zInteractionPoint,
                                           RandomNumberGenerator& rng) const = 0;
};

}  // namespace interactions
}  // namespace simprop

#endif
