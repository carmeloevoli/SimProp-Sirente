#include <limits>

#include "simprop/analyticalSolutions/beniamino.h"

namespace simprop {
namespace solutions {

#define INTSTEPS 1000

double Beniamino::I_dEpsilon(double EnuObserved, double Ep, double z) const {
  const auto x = EnuObserved * (1. + z) / Ep;
  const auto eta = 4. * Ep / pow2(SI::protonMassC2);
  const auto meanBkgPhotonEnergy = SI::kBoltzmann * 2.7 * SI::K * (1. + z);
  auto integrand = [this, Ep, eta, x, z](double lnEpsilon) {
    auto epsilon = std::exp(lnEpsilon);
    auto value = epsilon * m_cmb->density(epsilon / (1. + z));
    value *= sigma_nus.get(epsilon * eta, x);
    return value;
  };
  return utils::QAGIntegration<double>(integrand, std::log(0.01 * meanBkgPhotonEnergy),
                                       std::log(100. * meanBkgPhotonEnergy), INTSTEPS, 1e-4);
}

double Beniamino::I_dEp(double EnuObserved, double minEp, double maxEp, double z) const {
  const auto zMax = 3.;
  auto integrand = [this, EnuObserved, z, zMax](double lnEp) {
    auto Ep = std::exp(lnEp);
    auto value = computeFlux(Ep, z, zMax);
    value *= I_dEpsilon(EnuObserved, Ep, z);
    return value;
  };
  return utils::QAGIntegration<double>(integrand, std::log(minEp), std::log(maxEp), INTSTEPS, 1e-3);
}

// double Beniamino::I_z(double Enu, double zMax, double relError) const {
//   auto integrand = [this, Enu](double z) {
//     auto value = m_cosmology->dtdz(z);
//     value *= pow3(1. + z);
//     value *= I_Ep(Enu * (1. + z), 1e4 * Enu * (1. + z), z);
//     return value;
//   };
//   return utils::QAGIntegration<double>(integrand, 0., zMax, INTSTEPS, relError);
// }

double Beniamino::computeNeutrinoFlux(double Enu, double zMax, double relError) const {
  auto value = SI::cLight / 4. / M_PI;
  // value *= I_z(Enu, zMax, relError);
  return value;
}

}  // namespace solutions
}  // namespace simprop