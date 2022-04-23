#include "simprop/energyLosses/BGG2002ContinuousLosses.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

BGG2002ContinuousLosses::BGG2002ContinuousLosses() : ContinuousLosses() {
  LOGD << "calling " << __func__ << " constructor";
}

double BGG2002ContinuousLosses::getInterpolated(double E) const {
  double b_l = 0;
  const auto logE = std::log10(E / SI::eV);
  if (m_totalLosses.xIsInside(logE)) {
    b_l = std::pow(10., m_totalLosses.get(logE));
    b_l /= SI::year;
  }
  return b_l;
}

double BGG2002ContinuousLosses::dlnGamma_dt(PID pid, double Gamma, double z) const {
  const auto E = Gamma * SI::protonMassC2 * getNucleusMassNumber(pid);
  const auto redshiftedEnergy = E * (1. + z);
  double b_l = getInterpolated(redshiftedEnergy);
  if (b_l > 0.) {
    b_l *= pow3(1. + z);
    const double Z = (double)getNucleusCharge(pid);
    const double A = (double)getNucleusCharge(pid);
    b_l *= pow2(Z) / A;
  }
  return std::max(b_l, 0.);
}

// double BGG2002ContinuousLosses::evolve(double E_i, double z_i, double z_f, PID pid) const {
//   if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
//   const auto b = dlnGamma_dz(E_i, z_i, pid);
//   if (b < 0) throw std::runtime_error("b must be positive");
//   const auto factor = 1. - (z_i - z_f) * b;
//   if (factor > 1.) throw std::runtime_error("dz too large for continuous energy losses");
//   return E_i * factor;
// }

}  // namespace losses
}  // namespace simprop
