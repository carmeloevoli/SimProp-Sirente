#ifndef SIMPROP_INTERACTIONS_PPPCMB_H
#define SIMPROP_INTERACTIONS_PPPCMB_H

#include "simprop/interactions/Interaction.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace interactions {

class pppCmb final : public Interaction {
 protected:
  std::shared_ptr<photonfields::PhotonField> m_cmb;

 public:
  pppCmb(const std::shared_ptr<xsecs::CrossSection>& sigma)
      : Interaction(sigma), m_cmb(std::make_shared<photonfields::CMB>()) {}
  virtual ~pppCmb() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;
};

}  // namespace interactions
}  // namespace simprop

#endif
