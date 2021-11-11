#include "simprop/cosmology/cosmology.h"

#include <cmath>
#include <stdexcept>

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

}  // namespace cosmo
}  // namespace simprop