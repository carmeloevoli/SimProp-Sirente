#ifndef SIMPROP_LOSSES_ADIABATIC_CONTINUOUS_H
#define SIMPROP_LOSSES_ADIABATIC_CONTINUOUS_H

#include "simprop/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/units.h"

namespace simprop {
namespace losses {

class AdiabaticContinuousLosses final : public ContinuousLosses {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;

 public:
  AdiabaticContinuousLosses(const std::shared_ptr<cosmo::Cosmology>& cosmology);
  virtual ~AdiabaticContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;
};

}  // namespace losses
}  // namespace simprop

#endif