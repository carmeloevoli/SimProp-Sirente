#include "simprop/particleStacks/SourceEvolutionBuilder.h"

#include "simprop/common.h"

namespace simprop {

SourceEvolutionBuilder::SourceEvolutionBuilder(PID pid, SourceEvolutionParams params, size_t size)
    : Builder(pid, size) {
  m_GammaRange = params.GammaRange;
  m_zRange = params.zRange;
  m_slope = params.slope;
  m_GammaCutoff = params.GammaCutoff;
  m_evolutionIndex = params.evolutionIndex;
  auto GammaMin = m_GammaRange.first;
  m_maxWeight = std::pow(GammaMin, -m_slope + 1.) * std::exp(-GammaMin / m_GammaCutoff);
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack SourceEvolutionBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto z_i = getRndRedshift(m_zRange, RandomNumber(rng()));
    auto Gamma_i = getRndGamma(m_GammaRange, RandomNumber(rng()));
    auto w_i = std::pow(Gamma_i, -m_slope + 1.) * std::exp(-Gamma_i / m_GammaCutoff);
    w_i *= std::pow(1. + z_i, m_evolutionIndex + 2.);
    stack.emplace_back(Particle{m_pid, Redshift(z_i), LorentzFactor(Gamma_i), w_i / m_maxWeight});
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