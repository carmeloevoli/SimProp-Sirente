#ifndef SIMPROP_CONTINUOUS_LOSSES_H
#define SIMPROP_CONTINUOUS_LOSSES_H

#include <memory>

#include "simprop/cosmology.h"
#include "simprop/pid.h"

namespace simprop {
namespace losses {

class ContinuousLosses {
 public:
  explicit ContinuousLosses(const std::shared_ptr<cosmo::Cosmology>& cosmology)
      : m_cosmology(cosmology) {}
  virtual ~ContinuousLosses() = default;

  virtual double dlnGamma_dz(PID pid, double Gamma, double z = 0) const = 0;
  virtual double dlnGamma_dt(PID pid, double Gamma, double z = 0) const = 0;

 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
};

}  // namespace losses
}  // namespace simprop

#endif