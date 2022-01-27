#ifndef SIMPROP_XSECS_ABSTRACTCROSSSECTION_H
#define SIMPROP_XSECS_ABSTRACTCROSSSECTION_H

#include "simprop/pid.h"

namespace simprop {
namespace xsecs {

class CrossSection {
 public:
  CrossSection() {}
  virtual ~CrossSection() = default;
  virtual double getAtEpsPrime(PID pid, double epsPrime) const = 0;
  virtual double getAtS(PID pid, double s) const = 0;
  virtual double getPhotonEnergyThreshold() const = 0;
};

}  // namespace xsecs
}  // namespace simprop

#endif