#ifndef SIMPROP_XSECS_ABSTRACTCROSSSECTION_H
#define SIMPROP_XSECS_ABSTRACTCROSSSECTION_H

#include "simprop/pid.h"

namespace simprop {
namespace xsecs {

class CrossSection {
 public:
  CrossSection() {}
  virtual ~CrossSection() = default;
  virtual double get(PID pid, double photonEnergy) const = 0;
  virtual double getThreshold() const = 0;
};

}  // namespace xsecs
}  // namespace simprop

#endif