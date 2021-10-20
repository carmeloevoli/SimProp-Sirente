#include "simprop/common.h"

#include <cassert>
#include <cmath>

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double r) {
  assert(!(r < 0. || r > 1.));
  assert(energyRange.second > energyRange.first);

  const auto C = std::exp(r * std::log(energyRange.second / energyRange.first));
  return C * energyRange.first;
}

double GetRndRedshift(double maxRedshift, int evolutionIndex, double r) {
  assert(!(r < 0. || r > 1.));
  assert(maxRedshift > 0.);
  assert(evolutionIndex > 1);

  const double n = (double)evolutionIndex;
  const auto C = r * (std::pow(1. + maxRedshift, 1. - n) - 1.);
  return std::pow(C + 1., 1. / (1. - n)) - 1.;
}

}  // namespace simprop