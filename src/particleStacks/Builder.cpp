#include "simprop/particleStacks/Builder.h"

namespace simprop {

Range Builder::getRedshiftRange(const ParticleStack& stack) const {
  auto r = std::minmax_element(
      stack.begin(), stack.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  return {r.first->getRedshift(), r.second->getRedshift()};
}

Range Builder::getGammaRange(const ParticleStack& stack) const {
  auto r = std::minmax_element(
      stack.begin(), stack.end(),
      [](const Particle& a, const Particle& b) { return a.getGamma() < b.getGamma(); });
  return {r.first->getGamma(), r.second->getGamma()};
}

}  // namespace simprop