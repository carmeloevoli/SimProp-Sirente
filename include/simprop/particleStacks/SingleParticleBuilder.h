#ifndef SIMPROP_SINGLEPARTICLEBUILDER_H
#define SIMPROP_SINGLEPARTICLEBUILDER_H

#include "simprop/particleStacks/Builder.h"

namespace simprop {

struct SingleParticleParams {
  double Gamma;
  double z;
};

class SingleParticleBuilder final : public Builder {
 protected:
  double m_Gamma = 1e12;
  double m_z = 1;

 public:
  SingleParticleBuilder(PID pid, SingleParticleParams params, size_t size = 1);
  ParticleStack build(RandomNumberGenerator& rng) const override;
  ParticleStack build() const;
};

}  // namespace simprop

#endif