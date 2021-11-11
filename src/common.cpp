#include "simprop/common.h"

#include <cassert>
#include <cmath>

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double slope, double r) {
  assert(!(r < 0. || r > 1.));
  assert(energyRange.second > energyRange.first);
  assert(slope >= 1.0);
  using std::exp;
  using std::log;

  if (slope == 1.0) {
    const auto C = exp(r * log(energyRange.second / energyRange.first));
    return C * energyRange.first;
  } else {
    const auto ER = 1. - r * (1. - pow(energyRange.second / energyRange.first, 1. - slope));
    const auto C = pow(ER, 1. / (1. - slope));
    return C * energyRange.first;
  }
}

double GetRndRedshift(std::pair<double, double> redshiftRange, int evolutionIndex, double r) {
  assert(!(r < 0. || r > 1.));
  assert(evolutionIndex > 1);
  using std::pow;

  const double n = (double)evolutionIndex;
  const auto ZM = pow(1. + redshiftRange.second, 1. - n);
  const auto Zm = pow(1. + redshiftRange.first, 1. - n);
  const auto C = pow(r * (ZM - Zm) + Zm, 1. / (1. - n));
  return C - 1;
}

}  // namespace simprop