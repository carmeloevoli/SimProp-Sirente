#ifndef SIMPROP_MONOCHROMATICBUILDER_H
#define SIMPROP_MONOCHROMATICBUILDER_H

#include "simprop/common.h"
#include "simprop/particleStacks/Builder.h"

namespace simprop {

struct MonochromaticParams {
  double Gamma;
  Range zRange;
  double evolutionIndex;
};

class MonochromaticBuilder final : public Builder {
 protected:
  double m_Gamma = -1;
  Range m_zRange = {0., 1.};
  double m_evolutionIndex = 1;
  double m_maxWeight;

 public:
  MonochromaticBuilder(PID pid, MonochromaticParams params, size_t size = 1);
  ParticleStack build(RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif