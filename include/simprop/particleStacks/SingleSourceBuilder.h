#ifndef SIMPROP_SINGLESOURCEBUILDER_H
#define SIMPROP_SINGLESOURCEBUILDER_H

#include "simprop/particleStacks/Builder.h"

namespace simprop {

class SingleSourceBuilder final : public Builder {
 protected:
  double m_z = 1;
  Range m_GammaRange = {1e8, 1e14};
  double m_slope = 2.;

 public:
  SingleSourceBuilder(PID pid, size_t size = 1);
  void setRedshift(double z) { m_z = z; };
  void setGammaRange(Range GammaRange) { m_GammaRange = GammaRange; };
  void setSlope(double slope) { m_slope = slope; }

  ParticleStack build(RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif