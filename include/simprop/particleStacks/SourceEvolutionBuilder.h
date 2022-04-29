#ifndef SIMPROP_SOURCEEVOLUTIONBUILDER_H
#define SIMPROP_SOURCEEVOLUTIONBUILDER_H

#include "simprop/common.h"
#include "simprop/particleStacks/Builder.h"

namespace simprop {

struct SourceEvolutionParams {
  Range GammaRange;
  Range zRange;
  double slope;
  double GammaCutoff;
  double evolutionIndex;
};

class SourceEvolutionBuilder final : public Builder {
 protected:
  Range m_GammaRange = {1e8, 1e14};
  Range m_zRange = {0., 1.};
  double m_slope = 2;
  double m_GammaCutoff = 1e14;
  double m_evolutionIndex = 1;
  double m_maxWeight;

 public:
  SourceEvolutionBuilder(PID pid, SourceEvolutionParams params, size_t size = 1);
  ParticleStack build(RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif