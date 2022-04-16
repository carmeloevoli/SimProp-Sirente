#include "simprop/particleStack.h"

#include "simprop/common.h"
#include "simprop/utils/logging.h"

namespace simprop {

double getRndGamma(Range gammaRange, double slope, double r) {
  using std::exp;
  using std::log;
  using std::pow;

  if (slope == 1.0) {
    const auto C = exp(r * log(gammaRange.second / gammaRange.first));
    return C * gammaRange.first;
  } else {
    const auto ER = 1. - r * (1. - pow(gammaRange.second / gammaRange.first, 1. - slope));
    const auto C = pow(ER, 1. / (1. - slope));
    return C * gammaRange.first;
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

ParticleStack::ParticleStack(PID pid, int nParticles, RandomNumberGenerator& rng)
    : m_pid(pid), m_size(nParticles), m_rng(rng) {}

void ParticleStack::buildSingleParticleStack(double z, double Gamma) {
  for (size_t i = 0; i < m_size; ++i) m_particles.push_back(Particle{m_pid, z, Gamma});
  assert(m_particles.size() == m_size);
}

void ParticleStack::buildMultipleParticleStack(Range zRange, Range gammaRange, double slope) {
  for (size_t i = 0; i < m_size; ++i) {
    auto z_i = getRndRedshift(zRange, 2, m_rng());
    auto Gamma_i = getRndGamma(gammaRange, slope, m_rng());
    m_particles.push_back(Particle{m_pid, z_i, Gamma_i});
  }
  LOGD << "built primaries with size " << m_particles.size();
  auto z_r = getRedshiftRange();
  LOGD << "z range (" << z_r.first << "," << z_r.second << ")";
  auto G_r = getGammaRange();
  LOGD << "Gamma range (" << G_r.first << "," << G_r.second << ")";
}

void ParticleStack::buildInitialState(Range zRange, Range gammaRange, double slope) {
  if (m_size == 1)
    buildSingleParticleStack(zRange.first, gammaRange.first);
  else
    buildMultipleParticleStack(zRange, gammaRange, slope);
}

Range ParticleStack::getRedshiftRange() const {
  auto r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  return {r.first->getRedshift(), r.second->getRedshift()};
}

Range ParticleStack::getGammaRange() const {
  auto r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getGamma() < b.getGamma(); });
  return {r.first->getGamma(), r.second->getGamma()};
}

}  // namespace simprop