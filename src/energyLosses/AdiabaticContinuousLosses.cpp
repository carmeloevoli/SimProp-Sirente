#include "simprop/energyLosses/AdiabaticContinuousLosses.h"

namespace simprop {
namespace losses {

double AdiabaticContinuousLosses::dlnE_dz(PID pid, double E, double z) const {
  const auto b_a = 1. / (1. + z);
  return b_a;
}

double AdiabaticContinuousLosses::dlnE_dt(PID pid, double E, double z) const {
  const auto b_a = 1. / (1. + z);
  return b_a / m_cosmology.dtdz(z);
}

// double AdiabaticContinuousLosses::evolve(double E_i, double z_i, double z_f, PID pid) const {
//   if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
//   return E_i * (1. + z_f) / (1. + z_i);
// }

}  // namespace losses
}  // namespace simprop
