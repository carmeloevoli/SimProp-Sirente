// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_CROSSSECTIONS_CROSSSECTION_H_
#define SIMPROP_CROSSSECTIONS_CROSSSECTION_H_

#include "simprop/core/pid.h"

namespace simprop {
namespace xsecs {

class CrossSection {
 public:
  CrossSection() {}
  virtual ~CrossSection() = default;
  virtual double getAtEpsPrime(PID pid, double eps) const = 0;
  virtual double getEpsPrimeThreshold() const = 0;
};

}  // namespace xsecs
}  // namespace simprop

#endif  // SIMPROP_CROSSSECTIONS_CROSSSECTION_H_
