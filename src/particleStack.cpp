#include "simprop/particleStack.h"

#include "simprop/common.h"
#include "simprop/utils/logging.h"

namespace simprop {

double getRndEnergy(std::pair<double, double> energyRange, double slope, double r) {
  using std::exp;
  using std::log;
  using std::pow;

  if (slope == 1.0) {
    const auto C = exp(r * log(energyRange.second / energyRange.first));
    return C * energyRange.first;
  } else {
    const auto ER = 1. - r * (1. - pow(energyRange.second / energyRange.first, 1. - slope));
    const auto C = pow(ER, 1. / (1. - slope));
    return C * energyRange.first;
  }
}

double getRndRedshift(std::pair<double, double> redshiftRange, int evolutionIndex, double r) {
  using std::pow;

  const double n = (double)evolutionIndex;
  const auto ZM = pow(1. + redshiftRange.second, 1. - n);
  const auto Zm = pow(1. + redshiftRange.first, 1. - n);
  const auto C = pow(r * (ZM - Zm) + Zm, 1. / (1. - n));
  return C - 1;
}

ParticleStack::ParticleStack(PID pid, int nParticles, int seed)
    : m_pid(pid), m_size(nParticles), m_rng(RandomNumberGenerator(seed)) {}

void ParticleStack::buildInitialStates(Range zRange, Range eRange, double slope) {
  for (size_t i = 0; i < m_size; ++i) {
    auto z_i = getRndRedshift(zRange, 2, m_rng());
    auto E_i = getRndEnergy(eRange, slope, m_rng());
    m_particles.push_back(Particle{m_pid, z_i, E_i});
  }
  LOGD << "built primaries with size " << m_particles.size();
  auto z_r = getRedshiftRange();
  LOGD << "z range (" << z_r.first << "," << z_r.second << ")";
  auto E_r = getEnergyRange();
  LOGD << "E range (" << E_r.first / SI::eV << "," << E_r.second / SI::eV << ")";
}

Range ParticleStack::getRedshiftRange() const {
  auto r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  return {r.first->getRedshift(), r.second->getRedshift()};
}

Range ParticleStack::getEnergyRange() const {
  auto r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getEnergy() < b.getEnergy(); });
  return {r.first->getEnergy(), r.second->getEnergy()};
}

}  // namespace simprop