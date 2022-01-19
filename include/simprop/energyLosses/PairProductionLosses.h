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

class PairProductionLosses final : public ContinuousLosses {
 protected:
  std::vector<std::shared_ptr<photonfields::PhotonField>> m_photonFields;

 public:
  explicit PairProductionLosses(
      const std::shared_ptr<cosmo::Cosmology>& cosmology,
      const std::vector<std::shared_ptr<photonfields::PhotonField>>& photonFields)
      : ContinuousLosses(cosmology), m_photonFields(photonFields) {}
  virtual ~PairProductionLosses() = default;

  double dlnGamma_dt(PID pid, double Gamma, double z = 0) const override;
  double dlnGamma_dz(PID pid, double Gamma, double z = 0) const override;
};

}  // namespace losses

}  // namespace simprop

#endif