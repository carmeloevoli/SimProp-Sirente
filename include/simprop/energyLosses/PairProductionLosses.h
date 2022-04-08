#ifndef SIMPROP_LOSSES_PAIRPRODUCTION_H
#define SIMPROP_LOSSES_PAIRPRODUCTION_H

#include <memory>
#include <vector>

#include "simprop/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/PhotonField.h"
#include "simprop/units.h"

namespace simprop {
namespace losses {

using PhotonFields = std::vector<std::shared_ptr<photonfields::PhotonField>>;

class PairProductionLosses final : public ContinuousLosses {
 protected:
  PhotonFields m_photonFields;

 public:
  explicit PairProductionLosses(const PhotonFields& photonFields)
      : ContinuousLosses(), m_photonFields(photonFields) {}
  virtual ~PairProductionLosses() = default;
  double dlnGamma_dt(PID pid, double Gamma, double z = 0) const override;

 protected:
  double dotGamma(double Gamma) const;
};

}  // namespace losses
}  // namespace simprop

#endif