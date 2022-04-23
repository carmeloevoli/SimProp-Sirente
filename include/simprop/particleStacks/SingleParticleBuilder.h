#ifndef SIMPROP_SINGLEPARTICLEBUILDER_H
#define SIMPROP_SINGLEPARTICLEBUILDER_H

#include "simprop/particleStacks/Builder.h"

namespace simprop {

class SingleParticleBuilder final : public Builder {
 protected:
  double m_z = 1;
  double m_Gamma = 1e12;

 public:
  SingleParticleBuilder(PID pid, size_t size = 1);
  void setRedshift(double z) { m_z = z; };
  void setGamma(double Gamma) { m_Gamma = Gamma; };

  ParticleStack build(RandomNumberGenerator& rng) const override;
  ParticleStack build() const;
};

}  // namespace simprop

#endif