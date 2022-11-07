#ifndef SIMPROP_LOSSES_BGG2006_CONTINUOUS_H
#define SIMPROP_LOSSES_BGG2006_CONTINUOUS_H

#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace losses {

class BGG2006ContinuousLosses final : public ContinuousLosses {
 protected:
  // Berezinsky, Gazizov & Grigorieva, 2006, PRD, vol. 74, Issue 4, id. 043005
  const std::string totalLossesFilename = "data/losses_pair_BGG2006.txt";
  utils::LookupArray<501> m_totalLosses;

 public:
  BGG2006ContinuousLosses();
  virtual ~BGG2006ContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double getInterpolated(double E) const;
};

}  // namespace losses
}  // namespace simprop

#endif