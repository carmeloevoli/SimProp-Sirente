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
  return H0 * sqrt(OmegaM * pow(x, 3.) + OmegaR * pow(x, 4.) + OmegaL);
}

double t_H(double z) { return 1.0 / H(z); }

double dtdz(double z) {
  assert(!(z < 0));
  using std::pow;
  using std::sqrt;

  const auto x = sqrt(OmegaL / OmegaM) * pow(1 + z, -3.0 / 2.0);
  const auto dxdz = sqrt(OmegaL / OmegaM) * pow(1 + z, -5.0 / 2.0) * (-3.0 / 2.0);
  const auto const1 = 2 * sqrt(1 + OmegaM / OmegaL) / (3.0 * H0);

  const auto numer = dxdz * (1 + x * pow(pow(x, 2) + 1, -0.5));
  const auto denom = x + sqrt(pow(x, 2) + 1);

  return (const1 * numer / denom);
}

double adiabaticLossesEnergyFraction(double z_i, double z_f) {
  assert(z_f < z_i);
  auto result =
      gsl::QAGIntegration<double>([](double z) { return -H(z) * dtdz(z); }, z_i, z_f, LIMIT);
  return std::exp(result);
}

}  // namespace cosmo
}  // namespace simprop