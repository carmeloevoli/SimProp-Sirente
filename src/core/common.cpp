#include "simprop/core/common.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace simprop {

double getRndLogUniform(Range range, double r) {
  assert(range.second >= range.first);
  if (range.second == range.first) {
    return range.second;
  } else {
    const auto C = std::pow(range.second / range.first, r);
    return C * range.first;
  }
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

double mu2t(double mu, double s) {
  constexpr auto mp2 = SI::protonMassC2 * SI::protonMassC2;
  constexpr auto mpi2 = SI::pionMassC2 * SI::pionMassC2;
  const auto sqrts = std::sqrt(s);
  const auto A = mp2 - (s + mp2) * (s + mp2 - mpi2) / 2. / s;
  const auto B = -(s - mp2) / sqrts * std::sqrt(pow2((s + mp2 - mpi2) / 2. / sqrts) - mp2);
  return A + B * mu;
}

double t2mu(double t, double s) {
  constexpr auto mp2 = SI::protonMassC2 * SI::protonMassC2;
  constexpr auto mpi2 = SI::pionMassC2 * SI::pionMassC2;
  const auto sqrts = std::sqrt(s);
  const auto A = mp2 - (s + mp2) * (s + mp2 - mpi2) / 2. / s;
  const auto B = -(s - mp2) / sqrts * std::sqrt(pow2((s + mp2 - mpi2) / 2. / sqrts) - mp2);
  return (t - A) / B;
}

}  // namespace simprop