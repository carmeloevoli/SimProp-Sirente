#include "simprop/particleStacks/SourceEvolutionBuilder.h"

namespace simprop {

ParticleStack SourceEvolutionBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto z_i = getRndRedshift(m_zRange, 2, rng());
    auto Gamma_i = getRndGamma(m_GammaRange, m_slope, rng());
    stack.emplace_back(Particle{m_pid, z_i, Gamma_i});
  }
  assert(stack.size() == m_size);
  LOGD << "built primaries with size " << stack.size();
  auto zRange = getRedshiftRange(stack);
  LOGD << "redshift in (" << zRange.first << "," << zRange.second << ")";
  auto GammaRange = getGammaRange(stack);
  LOGD << "Gamma in (" << GammaRange.first << "," << GammaRange.second << ")";
  return stack;
}

}  // namespace simprop