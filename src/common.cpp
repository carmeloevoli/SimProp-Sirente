#include "simprop/common.h"

#include <cassert>
#include <cmath>

namespace simprop {

double getRndGamma(Range gammaRange, RandomNumber r) {
  using std::pow;
  const auto C = pow(gammaRange.second / gammaRange.first, r.get());
  return C * gammaRange.first;
}

double getRndRedshift(Range redshiftRange, RandomNumber r) {
  using std::pow;
  const double p = 2.0;
  const auto Zm = pow(1. + redshiftRange.first, -p + 1.);
  const auto ZM = pow(1. + redshiftRange.second, -p + 1.);
  const auto C = pow(r.get() * ZM + (1. - r.get()) * Zm, 1. / (-p + 1.));
  return C - 1;
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