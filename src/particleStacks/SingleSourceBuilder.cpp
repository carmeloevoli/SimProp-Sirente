#include "simprop/particleStacks/SingleSourceBuilder.h"

#include "simprop/core/common.h"

namespace simprop {

SingleSourceBuilder::SingleSourceBuilder(PID pid, SingleSourceParams params, size_t size)
    : Builder(pid, size) {
  m_GammaRange = params.GammaRange;
  m_z = params.z;
  m_slope = params.slope;
  m_GammaCutoff = params.GammaCutoff;
  auto GammaMin = m_GammaRange.first;
  m_maxWeight = std::pow(GammaMin, -m_slope + 1.) * std::exp(-GammaMin / m_GammaCutoff);
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack SingleSourceBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto Gamma_i = getRndLogUniform(m_GammaRange, rng());
    auto w_i = std::pow(Gamma_i, -m_slope + 1.) * std::exp(-Gamma_i / m_GammaCutoff);
    stack.emplace_back(Particle{m_pid, m_z, Gamma_i, w_i / m_maxWeight});
  }
  assert(stack.size() == m_size);
  LOGD << "built primaries with size " << stack.size();
  auto GammaRange = getGammaRange(stack);
  LOGD << "Gamma in (" << GammaRange.first << "," << GammaRange.second << ")";
  return stack;
}

}  // namespace simprop