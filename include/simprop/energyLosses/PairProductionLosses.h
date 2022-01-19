#ifndef SIMPROP_LOSSES_PAIRPRODUCTION_H
#define SIMPROP_LOSSES_PAIRPRODUCTION_H

#include "simprop/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/units.h"

namespace simprop {
namespace losses {

class PairProductionLosses final : public ContinuousLosses {
 public:
  explicit PairProductionLosses(const cosmo::Cosmology& cosmology) : ContinuousLosses(cosmology) {}
  virtual ~PairProductionLosses() = default;

  double dlnE_dt(PID pid, double E, double z = 0) const override;
  double dlnE_dz(PID pid, double E, double z = 0) const override;
  // double evolve(double E_i, double z_i, double z_f, PID pid) const override;
  // double evolve_rk4(double E_i, double z_i, double z_f, PID pid) const;

 protected:
};

}  // namespace losses

}  // namespace simprop

#endif