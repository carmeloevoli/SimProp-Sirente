#include "simprop/particleStacks/SingleParticleBuilder.h"

namespace simprop {

ParticleStack SingleParticleBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  LOGD << "building stack with " << m_size << " particles";
  LOGD << "type = " << getPidName(m_pid) << ", z = " << m_z << ", Gamma = " << m_Gamma;
  for (size_t i = 0; i < m_size; ++i) stack.push_back(Particle{m_pid, m_z, m_Gamma, true});
  assert(stack.size() == m_size);
  return stack;
}

}  // namespace  simprop
