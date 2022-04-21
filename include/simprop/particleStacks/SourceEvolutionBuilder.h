#ifndef SIMPROP_SOURCEEVOLUTIONBUILDER_H
#define SIMPROP_SOURCEEVOLUTIONBUILDER_H

#include "simprop/particleStacks/Builder.h"

namespace simprop {

class SourceEvolutionBuilder final : public Builder {
 protected:
  Range m_zRange = {0., 1.};
  Range m_GammaRange = {1e8, 1e14};
  double m_slope = 2;
  double m_evolutionIndex = 1;

 public:
  SourceEvolutionBuilder(PID pid, size_t size = 1) : Builder(pid, size) {}
  void setRedshiftRange(Range zRange) { m_zRange = zRange; };
  void setGammaRange(Range GammaRange) { m_GammaRange = GammaRange; };
  void setSlope(double slope) { m_slope = slope; }
  void setEvolutionIndex(double index) { m_evolutionIndex = index; }

  ParticleStack build(RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif