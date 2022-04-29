#ifndef SIMPROP_SINGLESOURCEBUILDER_H
#define SIMPROP_SINGLESOURCEBUILDER_H

#include "simprop/common.h"
#include "simprop/particleStacks/Builder.h"

namespace simprop {

struct SingleSourceParams {
  Range GammaRange;
  double z;
  double slope;
  double GammaCutoff;
};

class SingleSourceBuilder final : public Builder {
 protected:
  Range m_GammaRange = {1e8, 1e14};
  double m_z = 1;
  double m_slope = -1;
  double m_GammaCutoff = -1;
  double m_maxWeight = -1;

 public:
  SingleSourceBuilder(PID pid, SingleSourceParams params, size_t size = 1);
  ParticleStack build(RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif