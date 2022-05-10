#ifndef SIMPROP_LOSSES_BGG2002_CONTINUOUS_H
#define SIMPROP_LOSSES_BGG2002_CONTINUOUS_H

#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace losses {

class BGG2002ContinuousLosses final : public ContinuousLosses {
 protected:
  const std::string totalLossesFilename = "data/losses_pair_BGG2002.txt";
  utils::LookupArray<501> m_totalLosses;

 public:
  BGG2002ContinuousLosses();
  virtual ~BGG2002ContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double getInterpolated(double E) const;
};

}  // namespace losses
}  // namespace simprop

#endif