#include "simprop/core/common.h"

#include <cassert>
#include <cmath>

namespace simprop {

double getRndLogUniform(Range range, double r) {
  assert(range.second > range.first);
  using std::pow;
  const auto C = pow(range.second / range.first, r);
  return C * range.first;
}

Range getRedshiftRange(const ParticleStack& stack) {
  auto r = std::minmax_element(
      stack.begin(), stack.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  return {r.first->getRedshift(), r.second->getRedshift()};
}

Range getGammaRange(const ParticleStack& stack) {
  auto r = std::minmax_element(
      stack.begin(), stack.end(),
      [](const Particle& a, const Particle& b) { return a.getGamma() < b.getGamma(); });
  return {r.first->getGamma(), r.second->getGamma()};
}

}  // namespace simprop