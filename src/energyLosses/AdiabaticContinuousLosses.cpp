#include "simprop/energyLosses/AdiabaticContinuousLosses.h"

#include "simprop/utils/io.h"

namespace simprop {
namespace losses {

AdiabaticContinuousLosses::AdiabaticContinuousLosses(
    const std::shared_ptr<cosmo::Cosmology>& cosmology)
    : ContinuousLosses(), m_cosmology(cosmology) {
  LOGD << "calling " << __func__ << " constructor";
}

double AdiabaticContinuousLosses::beta(PID pid, double Gamma, double z) const {
  const auto b_a = 1. / (1. + z);
  return b_a / m_cosmology->dtdz(z);
}

// double AdiabaticContinuousLosses::evolve(double E_i, double z_i, double z_f, PID pid) const {
//   if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
//   return E_i * (1. + z_f) / (1. + z_i);
// }

}  // namespace losses
}  // namespace simprop
