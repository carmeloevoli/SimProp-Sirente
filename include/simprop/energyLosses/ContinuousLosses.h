#ifndef SIMPROP_CONTINUOUS_LOSSES_H
#define SIMPROP_CONTINUOUS_LOSSES_H

#include <memory>

#include "simprop/cosmology.h"
#include "simprop/pid.h"

namespace simprop {
namespace losses {

class ContinuousLosses {
 public:
  ContinuousLosses() {}
  virtual ~ContinuousLosses() = default;
  virtual double beta(PID pid, double Gamma, double z = 0) const = 0;
};

}  // namespace losses
}  // namespace simprop

#endif