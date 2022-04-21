#include "simprop/particleStacks/SingleParticleBuilder.h"

namespace simprop {

ParticleStack SingleParticleBuilder::build(RandomNumberGenerator& rng) const {
  ParticleStack stack;
  stack.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) stack.emplace_back(Particle{m_pid, m_z, m_Gamma});
  assert(stack.size() == m_size);
  return stack;
}

}  // namespace  simprop
