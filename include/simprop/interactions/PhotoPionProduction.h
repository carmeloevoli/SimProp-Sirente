#ifndef SIMPROP_INTERACTIONS_PPPEBL_H
#define SIMPROP_INTERACTIONS_PPPEBL_H

#include <memory>

#include "simprop/core/units.h"
#include "simprop/crossSections/PhotoPionXsecs.h"
#include "simprop/interactions/Interaction.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

using ProtonXsec = std::shared_ptr<xsecs::PhotoPionProtonXsec>;
using NeutronXsec = std::shared_ptr<xsecs::PhotoPionNeutronXsec>;

double sampleS(double r, double sMax, const std::shared_ptr<xsecs::CrossSection>& xs);

double sampleEps(double r, double nucleonEnergy, double z,
                 const std::shared_ptr<xsecs::CrossSection>& xs,
                 const std::shared_ptr<photonfields::PhotonField>& phField);

class PhotoPionProduction final : public Interaction {
 protected:
  const double m_sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);
  std::pair<ProtonXsec, NeutronXsec> m_xsecs;
  // std::shared_ptr<xsecs::PhotoPionProtonXsec> m_xs_proton;
  // std::shared_ptr<xsecs::PhotoPionNeutronXsec> m_xs_neutron;

 public:
  PhotoPionProduction(const std::shared_ptr<photonfields::PhotonField>& phField);
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;

  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;
  // TODO make protected:
  double computeRateComoving(double Gamma, double z,
                             const std::shared_ptr<xsecs::CrossSection>& xs) const;
};

}  // namespace interactions
}  // namespace simprop

#endif