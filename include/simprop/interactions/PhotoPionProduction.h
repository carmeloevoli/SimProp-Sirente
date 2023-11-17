#ifndef SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTION_H
#define SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTION_H

#include <fstream>
#include <iostream>
#include <memory>

#include "simprop/core/units.h"
#include "simprop/crossSections/PhotoPionXsecs.h"
#include "simprop/interactions/Interaction.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

PID pickNucleon(double r, PID pid);

class PhotoPionProduction : public Interaction {
 protected:
  const double m_sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  xsecs::PhotoPionXsec m_xs;

 public:
  PhotoPionProduction(const std::shared_ptr<photonfields::PhotonField>& phField);
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;
  double sampleS(double r, PID pid, double sMax) const;
  double sampleEps(double r, PID nucleon, double nucleonEnergy, double z) const;

  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;

 protected:
  double epsPdfIntegral(double photonEnergy, PID nucleon, double nucleonEnergy, double z) const;
};

}  // namespace interactions
}  // namespace simprop

#endif  // SIMPROP_INTERACTIONS_PHOTOPIONPRODUCTION_H
