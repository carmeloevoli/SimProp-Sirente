#ifndef SIMPROP_ENERGYLOSSES_ADIABATIC_CONTINUOUS_H
#define SIMPROP_ENERGYLOSSES_ADIABATIC_CONTINUOUS_H

#include "simprop/cosmology/cosmology.h"
#include "simprop/energyLosses/AbstractContinuousLosses.h"
#include "simprop/units.h"

namespace simprop {
namespace losses {

class AdiabaticContinuousLosses : public AbstractContinuousLosses {
 public:
  AdiabaticContinuousLosses() {}
  virtual ~AdiabaticContinuousLosses() = default;

  double dlnGamma_dz(double z, double E, PID pid) const override;
  double evolve(double E_i, double z_i, double z_f, PID pid) const override;
  // double evolve_rk4(double E_i, double z_i, double z_f, PID pid) const;
};

}  // namespace losses
}  // namespace simprop

#endif