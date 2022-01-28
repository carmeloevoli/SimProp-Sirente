#ifndef SIMPROP_INTERACTIONS_PPPEBL_H
#define SIMPROP_INTERACTIONS_PPPEBL_H

#include <memory>

#include "simprop/interactions/Interaction.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace interactions {

class PhotoPionProduction final : public Interaction {
 protected:
  std::shared_ptr<photonfields::PhotonField> m_ebl;

 public:
  PhotoPionProduction(const std::shared_ptr<xsecs::CrossSection>& sigma,
                      const std::shared_ptr<photonfields::PhotonField>& ebl)
      : Interaction(sigma), m_ebl(ebl) {}
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;

 protected:
  double computeRateComoving(PID pid, double Gamma, double z) const;
};

}  // namespace interactions
}  // namespace simprop

#endif