#ifndef SIMPROP_ENERGYLOSSES_ABSTRACTCONTINUOUS_H
#define SIMPROP_ENERGYLOSSES_ABSTRACTCONTINUOUS_H

#include "simprop/pid.h"

namespace simprop {
namespace losses {

class AbstractContinuousLosses {
 public:
  AbstractContinuousLosses() {}
  virtual ~AbstractContinuousLosses() = default;

  virtual double dlnGamma_dz(double z, double E, PID pid) const = 0;
  virtual double evolve(double E_i, double z_i, double z_f, PID pid) const = 0;
};

}  // namespace losses

}  // namespace simprop

#endif