#include "simprop/particleStacks/SingleSourceBuilder.h"

namespace simprop {

ParticleStack SingleSourceBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto Gamma_i = getRndGamma(m_GammaRange, m_slope, rng());
    stack.emplace_back(Particle{m_pid, m_z, Gamma_i, true});
  }
  assert(stack.size() == m_size);
  LOGD << "built primaries with size " << stack.size();
  auto GammaRange = getGammaRange(stack);
  LOGD << "Gamma in (" << GammaRange.first << "," << GammaRange.second << ")";
  return stack;
}

}  // namespace simprop