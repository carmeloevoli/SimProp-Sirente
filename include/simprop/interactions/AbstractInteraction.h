#ifndef SIMPROP_INTERACTIONS_ABSTRACTINTERACTION_H
#define SIMPROP_INTERACTIONS_ABSTRACTINTERACTION_H

#include "simprop/pid.h"

namespace simprop {
namespace interactions {

class AbstractInteration {
 public:
  AbstractInteration() {}
  virtual ~AbstractInteration() = default;

  virtual double getSigma(PID pid, double photonEnergy) const = 0;
};

}  // namespace interactions
}  // namespace simprop

#endif