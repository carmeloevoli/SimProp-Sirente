#ifndef SIMPROP_INTERACTIONS_PPPEBL_H
#define SIMPROP_INTERACTIONS_PPPEBL_H

#include <memory>

#include "simprop/core/units.h"
#include "simprop/crossSections/PhotoPionXsecs.h"
#include "simprop/interactions/Interaction.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

class PhotoPionProduction final : public Interaction {
 protected:
  const double m_sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  xsecs::PhotoPionXsec m_xs;

 public:
  PhotoPionProduction(const std::shared_ptr<photonfields::PhotonField>& phField);
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;

  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;

 public:  // TODO to be changed in private
  double computeRateComoving(PID pid, double Gamma, double z) const;
  double sampleS(double r, PID pid, double sMax) const;
  double epsPdfIntegral(double photonEnergy, PID nucleon, double nucleonEnergy, double z) const;
  double sampleEps(double r, PID nucleon, double nucleonEnergy, double z) const;
};

}  // namespace interactions
}  // namespace simprop

#endif