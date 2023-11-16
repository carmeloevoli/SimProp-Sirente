#include "simprop/crossSections/BreitWheeler.h"

#include <cmath>

#include "simprop/core/units.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace BreitWheeler {

double sigmaInCoMFrame(const double &s) {
  const auto chi = s / 4. / pow2(SI::electronMassC2);

  if (chi < 1. || chi > 1e5) return 0.;

  const auto beta = sqrt(1. - 1. / chi);
  return 3. / 16. * SI::sigmaTh * (1. - pow2(beta)) *
         (2. * beta * (pow2(beta) - 2.) + (3. - pow4(beta)) * log((1. + beta) / (1. - beta)));
}

double sigma(const double &eGamma, const double &eBkg, const double &mu) {
  using std::log;
  using std::sqrt;

  const auto s = 2. * eGamma * eBkg * (1. - mu);
  return sigmaInCoMFrame(s);
}

}  // namespace BreitWheeler
}  // namespace simprop
