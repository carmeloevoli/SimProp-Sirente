#ifndef SIMPROP_LOSSES_PAIRPRODUCTION_H
#define SIMPROP_LOSSES_PAIRPRODUCTION_H

#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/PhotonField.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace losses {

class PairProductionLosses final : public ContinuousLosses {
 protected:
  photonfields::PhotonFields m_photonFields;
  utils::LookupTable<2000, 100> m_betaProtons;
  bool m_doCaching = false;

 public:
  PairProductionLosses(const std::shared_ptr<photonfields::PhotonField>& photonField);
  PairProductionLosses(const photonfields::PhotonFields& photonFields);
  PairProductionLosses& doCaching();

  virtual ~PairProductionLosses() = default;
  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double computeProtonBeta(double Gamma, double z) const;
};

}  // namespace losses
}  // namespace simprop

#endif