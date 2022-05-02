#include "simprop/particleStacks/MonochromaticBuilder.h"

#include "simprop/common.h"

namespace simprop {

MonochromaticBuilder::MonochromaticBuilder(PID pid, MonochromaticParams params, size_t size)
    : Builder(pid, size) {
  m_Gamma = params.Gamma;
  m_zRange = params.zRange;
  m_evolutionIndex = params.evolutionIndex;
  m_maxWeight = 1.;
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack MonochromaticBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    auto z_i = getRndRedshift(m_zRange, RandomNumber(rng()));
    auto w_i = std::pow(1. + z_i, m_evolutionIndex + 1.);
    stack.emplace_back(Particle{m_pid, Redshift(z_i), LorentzFactor(m_Gamma), w_i / m_maxWeight});
  }
  assert(stack.size() == m_size);
  LOGD << "built primaries with size " << stack.size();
  auto zRange = getRedshiftRange(stack);
  LOGD << "redshift in (" << zRange.first << "," << zRange.second << ")";
  return stack;
}

}  // namespace simprop