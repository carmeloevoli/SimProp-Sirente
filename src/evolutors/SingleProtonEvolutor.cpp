#include "simprop/evolutors/singleProtonEvolutor.h"

#include <algorithm>

#include "simprop/utils/numeric.h"

#define VERYLARGEENERGY (1e25 * SI::eV)

namespace simprop {
namespace evolutors {

auto IsActive = [](const Particle& p) {
  constexpr double minPropagatingGamma = 1e6;
  return (p.isNucleus() && p.isActive() && p.getRedshift() > 1e-20 &&
          p.getGamma() > minPropagatingGamma);
};

SingleProtonEvolutor::SingleProtonEvolutor(RandomNumberGenerator& rng) : m_rng(rng) {}

double SingleProtonEvolutor::computeDeltaGamma(const Particle& particle,
                                               double deltaRedshift) const {
  const auto pid = particle.getPid();
  const auto Gamma = particle.getGamma();
  const auto zNow = particle.getRedshift();
  const auto zHalf = zNow - 0.5 * deltaRedshift;
  const auto zNext = zNow - deltaRedshift;
  double betaNow = 0, betaHalf = 0, betaNext = 0;
  for (auto losses : m_continuousLosses) {
    betaNow += losses->beta(pid, Gamma, zNow);
    betaHalf += losses->beta(pid, Gamma, zHalf);
    betaNext += losses->beta(pid, Gamma, zNext);
  }
  double value = deltaRedshift / 6.;
  value *= betaNow * m_cosmology->dtdz(zNow) + 4. * betaHalf * m_cosmology->dtdz(zHalf) +
           betaNext * m_cosmology->dtdz(zNext);
  return 1.0 - std::exp(-value);
}

double SingleProtonEvolutor::computeLossesRedshiftInterval(const Particle& particle) const {
  const auto zNow = particle.getRedshift();
  double dz = zNow;

  double deltaGamma = computeDeltaGamma(particle, zNow);
  if (deltaGamma > deltaGammaCritical) {
    dz = utils::rootFinder<double>(
        [&](double x) { return computeDeltaGamma(particle, x) - deltaGammaCritical; }, 0., zNow,
        100, 1e-3);
  }
  return dz;
}

void SingleProtonEvolutor::run(ParticleStack& stack) {
  size_t counter = 0;
  auto it = stack.begin();
  while (it != stack.end()) {
    const auto nowRedshift = it->getRedshift();
    const auto Gamma = it->getGamma();
    const auto dz_c = computeLossesRedshiftInterval(*it);
    assert(dz_c > 0. && dz_c <= nowRedshift);
    const auto deltaGamma = computeDeltaGamma(*it, dz_c);
    it->getNow() = {nowRedshift - dz_c, Gamma * (1. - deltaGamma)};
    it = std::find_if(it, stack.end(), IsActive);
    counter++;
  }
}  // run()

// double SingleProtonEvolutor::getObservedEnergy() const {
//   assert(m_stack.size() == 1);  // && m_stack[0].getRedshift() < 1e-20);
//   return m_stack[0].getGamma() * SI::protonMassC2;
// }

}  // namespace evolutors
}  // namespace simprop