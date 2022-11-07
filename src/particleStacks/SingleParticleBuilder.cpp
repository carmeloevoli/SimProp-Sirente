#include "simprop/particleStacks/SingleParticleBuilder.h"

namespace simprop {

SingleParticleBuilder::SingleParticleBuilder(PID pid, SingleParticleParams params, size_t size)
    : Builder(pid, size) {
  m_Gamma = params.Gamma;
  m_z = params.z;
  LOGD << "calling " << __func__ << " constructor";
}

ParticleStack SingleParticleBuilder::build() const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    stack.emplace_back(Particle{m_pid, m_z, m_Gamma});
  }
  assert(stack.size() == m_size);
  LOGD << "building stack with " << m_size << " particles";
  LOGD << "type = " << getPidName(m_pid) << ", z = " << m_z << ", Gamma = " << m_Gamma;
  return stack;
}

ParticleStack SingleParticleBuilder::build(RandomNumberGenerator& rng) const { return build(); }

}  // namespace  simprop
