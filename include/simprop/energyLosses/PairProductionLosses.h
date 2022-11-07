#ifndef SIMPROP_LOSSES_PAIRPRODUCTION_H
#define SIMPROP_LOSSES_PAIRPRODUCTION_H

#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace losses {

class PairProductionLosses final : public ContinuousLosses {
 protected:
  photonfields::PhotonFields m_photonFields;

 public:
  PairProductionLosses(const std::shared_ptr<photonfields::PhotonField>& photonField);
  PairProductionLosses(const photonfields::PhotonFields& photonFields);
  virtual ~PairProductionLosses() = default;
  double beta(PID pid, double Gamma, double z = 0) const override;

 public:
  double betaComoving(double Gamma) const;
};

}  // namespace losses
}  // namespace simprop

#endif