#include "simprop/particleStacks/SingleSourceBuilder.h"

#include "simprop/common.h"

namespace simprop {

SingleSourceBuilder::SingleSourceBuilder(PID pid, size_t size) : Builder(pid, size) {
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack SingleSourceBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto Gamma_i = getRndGamma(m_GammaRange, m_slope, rng());
    stack.emplace_back(Particle{m_pid, Redshift(m_z), LorentzFactor(Gamma_i), true});
  }
  assert(stack.size() == m_size);
  LOGD << "built primaries with size " << stack.size();
  auto GammaRange = getGammaRange(stack);
  LOGD << "Gamma in (" << GammaRange.first << "," << GammaRange.second << ")";
  return stack;
}

}  // namespace simprop