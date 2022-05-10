#ifndef SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H
#define SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H

#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace losses {

class PhotoPionContinuousLosses final : public ContinuousLosses {
 protected:
  utils::LookupArray<10000> m_totalLosses;
  const double m_inelasticity = 0.15;

 public:
  PhotoPionContinuousLosses();
  virtual ~PhotoPionContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double getInterpolated(double E) const;
};

}  // namespace losses

}  // namespace simprop

#endif