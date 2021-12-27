#ifndef SIMPROP_LOSSES_ADIABATIC_CONTINUOUS_H
#define SIMPROP_LOSSES_ADIABATIC_CONTINUOUS_H

#include "simprop/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/units.h"

namespace simprop {
namespace losses {

class AdiabaticContinuousLosses final : public ContinuousLosses {
 public:
  explicit AdiabaticContinuousLosses(const cosmo::Cosmology& cosmology)
      : ContinuousLosses(cosmology) {}
  virtual ~AdiabaticContinuousLosses() = default;

  double dlnE_dz(PID pid, double E, double z = 0) const override;
  double dlnE_dt(PID pid, double E, double z = 0) const override;
  // double evolve(double E_i, double z_i, double z_f, PID pid) const override;
  //  double evolve_rk4(double E_i, double z_i, double z_f, PID pid) const;
};

}  // namespace losses
}  // namespace simprop

#endif