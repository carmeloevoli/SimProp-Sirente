#include "simprop/evolutors/singleProtonEvolutor.h"

#include <algorithm>

#include "simprop/utils/logging.h"
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

double SingleProtonEvolutor::totalLosses(PID pid, double Gamma, double z) const {
  return std::accumulate(
      m_continuousLosses.begin(), m_continuousLosses.end(), 0.,
      [pid, Gamma, z](double beta, std::shared_ptr<losses::ContinuousLosses> losses) {
        return beta + losses->beta(pid, Gamma, z);
      });
}

double SingleProtonEvolutor::totalRate(PID pid, double Gamma, double z) const {
  return std::accumulate(
      m_interactions.begin(), m_interactions.end(), 0.,
      [pid, Gamma, z](double rate, std::shared_ptr<interactions::Interaction> interaction) {
        return rate + interaction->rate(pid, Gamma, z);
      });
}

double SingleProtonEvolutor::computeDeltaGamma(const Particle& particle,
                                               double deltaRedshift) const {
  const auto pid = particle.getPid();
  const auto Gamma = particle.getGamma();

  const auto zNow = particle.getRedshift();
  const auto zHalf = zNow - 0.5 * deltaRedshift;
  const auto zNext = zNow - deltaRedshift;

  const auto betaNow = totalLosses(pid, Gamma, zNow);
  const auto betaHalf = totalLosses(pid, Gamma, zHalf);
  const auto betaNext = totalLosses(pid, Gamma, zNext);

  // double betaNow = 0, betaHalf = 0, betaNext = 0;
  // for (auto losses : m_continuousLosses) {
  //   betaNow += losses->beta(pid, Gamma, zNow);
  //   betaHalf += losses->beta(pid, Gamma, zHalf);
  //   betaNext += losses->beta(pid, Gamma, zNext);
  // }

  double value = deltaRedshift / 6.;
  value *= betaNow * m_cosmology->dtdz(zNow) + 4. * betaHalf * m_cosmology->dtdz(zHalf) +
           betaNext * m_cosmology->dtdz(zNext);
  return 1.0 - std::exp(-value);
}

double SingleProtonEvolutor::computeLossesRedshiftInterval(const Particle& particle) const {
  const auto zNow = particle.getRedshift();
  const auto deltaGamma = computeDeltaGamma(particle, zNow);
  double dz = zNow;
  if (deltaGamma > deltaGammaCritical) {
    dz = utils::rootFinder<double>(
        [&](double x) { return computeDeltaGamma(particle, x) - deltaGammaCritical; }, 0., zNow,
        100, 1e-3);
  }
  return dz;
}

double SingleProtonEvolutor::computeInteractionRedshiftInterval(const Particle& particle) const {
  const auto pid = particle.getPid();
  const auto Gamma = particle.getGamma();
  const auto zNow = particle.getRedshift();
  const auto dt = std::fabs(1. / totalRate(pid, Gamma, zNow));
  // TODO why to put the fabs?
  return -dt / m_cosmology->dtdz(zNow) * std::log(1. - m_rng());
}

// void SingleProtonEvolutor::run(ParticleStack& stack) {
//   size_t counter = 0;
//   auto it = stack.begin();
//   while (it != stack.end()) {
//     const auto nowRedshift = it->getRedshift();
//     const auto Gamma = it->getGamma();
//     const auto dz_c = computeLossesRedshiftInterval(*it);
//     assert(dz_c > 0. && dz_c <= nowRedshift);
//     const auto deltaGamma = computeDeltaGamma(*it, dz_c);
//     it->getNow() = {nowRedshift - dz_c, Gamma * (1. - deltaGamma)};
//     it = std::find_if(it, stack.end(), IsActive);
//     counter++;
//   }
// }  // run()

void SingleProtonEvolutor::run(ParticleStack& stack) {
  size_t counter = 0;
  auto it = stack.begin();
  // size_t iniSize = stack.size();
  // const double initEnergy = sumEnergy();
  while (it != stack.end()) {
    const auto distance = it - stack.begin();
    const auto nowRedshift = it->getRedshift();
    const auto dz_s = computeInteractionRedshiftInterval(*it);
    const auto dz_c = computeLossesRedshiftInterval(*it);
    assert(dz_s > 0. && dz_c > 0. && dz_c <= nowRedshift);
    if (dz_s > dz_c || dz_s > nowRedshift) {
      const auto Gamma = it->getGamma();
      const auto dz = dz_c;
      const auto deltaGamma = computeDeltaGamma(*it, dz);
      it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
    } else {
      it->deactivate();
      const auto dz = dz_s;
      auto finalState = m_interactions[0]->finalState(*it, nowRedshift - dz, m_rng);
      stack.insert(stack.end(), finalState.begin(), finalState.end());
    }
    it = std::find_if(stack.begin() + distance, stack.end(), IsActive);
    counter++;
  }
}  // run()

// double SingleProtonEvolutor::getObservedEnergy() const {
//   assert(m_stack.size() == 1);  // && m_stack[0].getRedshift() < 1e-20);
//   return m_stack[0].getGamma() * SI::protonMassC2;
// }

}  // namespace evolutors
}  // namespace simprop