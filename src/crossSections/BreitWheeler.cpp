#include "simprop/crossSections/BreitWheeler.h"

#include <cmath>

#include "simprop/core/units.h"

namespace simprop {
namespace BreitWheeler {

double sigma(const double &eGamma, const double &eBkg, const double &mu) {
  using std::log;
  using std::sqrt;

  const auto s = 2. * eGamma * eBkg * (1. - mu);
  const auto chi = s / 4. / pow2(SI::electronMassC2);

  if (chi < 1.) return 0.;

  const auto beta = sqrt(1. - 1. / chi);
  return 3. / 16. * SI::sigmaTh * (1. - pow2(beta)) *
         (2. * beta * (pow2(beta) - 2.) + (3. - pow4(beta)) * log((1. + beta) / (1. - beta)));
}

}  // namespace BreitWheeler
}  // namespace simprop
