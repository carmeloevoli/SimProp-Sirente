#ifndef SIMPROP_INTERACTIONS_PPPEBL_H
#define SIMPROP_INTERACTIONS_PPPEBL_H

#include <memory>

#include "simprop/interactions/Interaction.h"
#include "simprop/photonFields/PhotonField.h"
#include "simprop/units.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

class PhotoPionProduction final : public Interaction {
 protected:
  std::shared_ptr<photonfields::PhotonField> m_ebl;
  const double m_sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);

 public:
  PhotoPionProduction(const std::shared_ptr<xsecs::CrossSection>& sigma,
                      const std::shared_ptr<photonfields::PhotonField>& ebl)
      : Interaction(sigma), m_ebl(ebl) {}
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;
  std::vector<Particle> finalState(PID pid, double Gamma, double z = 0) const override;

  double sample_s(double r, double sMax) const;
  double sample_eps(double r, double E, double z) const;
  double epsPdfIntegral(double photonEnergy, double E, double z) const;

  double computeRateComoving(PID pid, double Gamma, double z) const;  // TODO to make protected
  double phi(double s) const;

  double samplePionEnergy(double E, double z, RandomNumberGenerator& rng) const;
};

}  // namespace interactions
}  // namespace simprop

#endif