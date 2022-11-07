#ifndef SIMPROP_XSECS_ABSTRACTCROSSSECTION_H
#define SIMPROP_XSECS_ABSTRACTCROSSSECTION_H

#include "simprop/core/pid.h"

namespace simprop {
namespace xsecs {

class CrossSection {
 public:
  CrossSection() {}
  virtual ~CrossSection() = default;
  virtual double getAtEpsPrime(double epsPrime) const = 0;
  virtual double getAtS(double s) const = 0;
  virtual double getPhiAtS(double s) const = 0;
  virtual double getPhotonEnergyThreshold() const = 0;
};

}  // namespace xsecs
}  // namespace simprop

#endif