#include "simprop/cosmology/cosmology.h"

#include <cassert>
#include <cmath>

#include "simprop/cosmology/Planck2018.h"
#include "simprop/utils/gsl.h"

using namespace simprop::planck2018;

#define LIMIT 1000

namespace simprop {
namespace cosmo {

double H(double z) {
  using std::pow;
  using std::sqrt;

  const auto x = 1. + z;
  return H0 * sqrt(OmegaM * pow(x, 3.) + OmegaL);
}

double t_H(double z) { return 1.0 / H(z); }

double dtdz(double z) { return 1. / H(z) / (1. + z); }

double adiabaticRelativeLoss(double z_i, double z_f) {
  assert(z_f < z_i);
  auto result = std::log((z_f + 1.) / (z_i + 1.));
  return std::exp(result);
}

}  // namespace cosmo
}  // namespace simprop