#include "simprop/particleStack.h"

namespace simprop {

void ParticleStack::buildInitialStates() {
  m_particles.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = GetRndRedshift(m_params.redshiftRange, 2, m_rng());
    const auto E_i = GetRndEnergy(m_params.energyRange, 1, m_rng());
    auto p = Particle{m_params.pid, z_i, E_i};
    m_particles.emplace_back(p);
  }
  LOGD << "built primaries with size " << m_particles.size();
  printStateRanges();
}

void ParticleStack::printStateRanges() const {
  auto z_r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  LOGD << "z range (" << z_r.first->getRedshift() << "," << z_r.second->getRedshift() << ")";
  auto E_r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getEnergy() < b.getEnergy(); });
  using SI::eV;
  LOGD << "E range (" << E_r.first->getEnergy() / eV << "," << E_r.second->getEnergy() / eV << ")";
}
}  // namespace simprop