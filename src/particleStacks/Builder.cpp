#include "simprop/particleStacks/Builder.h"

namespace simprop {

double getRndGamma(Range gammaRange, double slope, double r) {
  using std::exp;
  using std::log;
  using std::pow;

  if (slope == 1.0) {
    const auto C = pow(gammaRange.second / gammaRange.first, r);
    return C * gammaRange.first;
  } else {
    const auto ER = 1. - r * (1. - pow(gammaRange.second / gammaRange.first, 1. - slope));
    const auto C = pow(ER, 1. / (1. - slope));
    return C * gammaRange.first;
  }
}

double getRndRedshift(Range redshiftRange, int evolutionIndex, double r) {
  using std::pow;

  const double m = (double)evolutionIndex;
  const auto ZM = pow(1. + redshiftRange.second, m + 1.);
  const auto Zm = pow(1. + redshiftRange.first, m + 1.);
  const auto C = pow(r * ZM + (1. - r) * Zm, 1. / (m + 1.));
  return C - 1;
}

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