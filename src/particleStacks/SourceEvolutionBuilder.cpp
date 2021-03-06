#include "simprop/particleStacks/SourceEvolutionBuilder.h"

#include "simprop/common.h"

namespace simprop {

SourceEvolutionBuilder::SourceEvolutionBuilder(PID pid, SourceEvolutionParams params,
                                               std::shared_ptr<cosmo::Cosmology> cosmology,
                                               size_t size)
    : Builder(pid, size), m_cosmology(cosmology) {
  m_GammaRange = params.GammaRange;
  m_zRange = params.zRange;
  m_slope = params.slope;
  // m_GammaCutoff = params.GammaCutoff;
  m_evolutionIndex = params.evolutionIndex;
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack SourceEvolutionBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  const auto y_min = (1. + m_zRange.first);
  const auto y_max = (1. + m_zRange.second);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = getRndLogUniform({y_min, y_max}, RandomNumber(rng())) - 1.;
    const auto Gamma_i = getRndLogUniform(m_GammaRange, RandomNumber(rng()));
    auto w_i = std::pow(Gamma_i / 1e8, 1. - m_slope);
    w_i *= std::pow(1. + z_i, m_evolutionIndex) / m_cosmology->E(z_i);
    stack.emplace_back(Particle{m_pid, Redshift(z_i), LorentzFactor(Gamma_i), w_i});
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