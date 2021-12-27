#ifndef SIMPROP_CONTINUOUS_LOSSES_H
#define SIMPROP_CONTINUOUS_LOSSES_H

#include <memory>

#include "simprop/cosmology.h"
#include "simprop/pid.h"

namespace simprop {
namespace losses {

class ContinuousLosses {
 public:
  explicit ContinuousLosses(const cosmo::Cosmology& cosmology) : m_cosmology(cosmology) {}
  virtual ~ContinuousLosses() = default;

  virtual double dlnE_dz(PID pid, double E, double z = 0) const = 0;
  virtual double dlnE_dt(PID pid, double E, double z = 0) const = 0;
  // virtual double evolve(double E_i, double z_i, double z_f, PID pid) const = 0;

 protected:
  const cosmo::Cosmology& m_cosmology;
};

}  // namespace losses
}  // namespace simprop

#endif