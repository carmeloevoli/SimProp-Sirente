#include "simprop/energyLosses/BGG2002ContinuousLosses.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

double BGG2002ContinuousLosses::dlnGamma_dt_0(double E, PID pid) const {
  double b_l = 0;
  const auto logE = std::log10(E / SI::eV);
  if (m_totalLosses.isWithinXRange(logE)) {
    b_l = std::pow(10., m_totalLosses.get(logE));
    b_l *= 1. / SI::year;
    const double Z = (double)getNucleusChargeNumber(pid);
    const double A = (double)getNucleusChargeNumber(pid);
    b_l *= utils::pow<2>(Z) / A;
  }
  return b_l;
}
double BGG2002ContinuousLosses::dlnGamma_dz(double z, double E, PID pid) const {
  const auto b_a = 1. / (1. + z);
  const auto redshiftedEnergy = E * (1. + z);
  double b_l = dlnGamma_dt_0(redshiftedEnergy, pid);
  if (b_l > 0.) {
    b_l *= utils::pow<3>(1. + z);
    b_l *= cosmo::dtdz(z);
  }
  return b_a + b_l;
}

double BGG2002ContinuousLosses::evolve(double E_i, double z_i, double z_f, PID pid) const {
  if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
  const auto b = dlnGamma_dz(E_i, z_i, pid);
  if (b < 0) throw std::runtime_error("b must be positive");
  const auto factor = 1. - (z_i - z_f) * b;
  if (factor > 1.) throw std::runtime_error("dz too large for continuous energy losses");
  return E_i * factor;
}

}  // namespace losses
}  // namespace simprop
