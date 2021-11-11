#include "simprop/energyLosses/continuousLosses.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace losses {

double ContinuousLosses::dlnGamma_dz(double z, double E, PID pid) const {
  const auto b_a = 1. / (1. + z);
  const auto redshiftedEnergy = E * (1. + z);
  double b_l = 0;
  if (m_totalLosses.isWithinXRange(redshiftedEnergy)) {
    auto b_l = std::pow(10., m_totalLosses.get(std::log10(redshiftedEnergy / SI::eV)));
    b_l *= 1. / SI::year;
    b_l *= utils::pow<3>(1. + z);
    b_l *= cosmo::dtdz(z);
    const double Z = (double)getNucleusChargeNumber(pid);
    const double A = (double)getNucleusChargeNumber(pid);
    b_l *= utils::pow<2>(Z) / A;
  }
  return b_a + b_l;
}

double ContinuousLosses::evolve(double E_i, double z_i, double z_f, PID pid) const {
  if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
  const auto b = dlnGamma_dz(E_i, z_i, pid);
  if (b < 0) throw std::runtime_error("b must be positive");
  const auto factor = 1. - (z_i - z_f) * b;
  if (factor > 1.) throw std::runtime_error("dz too large for continuous energy losses");
  return E_i * factor;
}

double ContinuousLosses::evolve_rk4(double E_i, double z_i, double z_f, PID pid) const {
  if (z_f > z_i) throw std::invalid_argument("z_f must be smaller than z_i");
  const auto h = z_i - z_f;
  const auto k_1 = E_i * dlnGamma_dz(E_i, z_i, pid);
  const auto k_2 = (E_i - .5 * h * k_1) * dlnGamma_dz(E_i - .5 * h * k_1, z_i - 0.5 * h, pid);
  const auto k_3 = (E_i - .5 * h * k_2) * dlnGamma_dz(E_i - .5 * h * k_2, z_i - 0.5 * h, pid);
  const auto k_4 = (E_i - h * k_3) * dlnGamma_dz(E_i - h * k_3, z_i - h, pid);
  return E_i - 1. / 6. * h * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
  // const auto b = dlnGamma_dz(E_i, z_i, pid);
  // if (b < 0) throw std::runtime_error("b must be positive");
  // const auto factor = 1. - (z_i - z_f) * b;
  // if (factor > 1.) throw std::runtime_error("dz too large for continuous energy losses");
  // return E_i * factor;
}

}  // namespace losses
}  // namespace simprop
