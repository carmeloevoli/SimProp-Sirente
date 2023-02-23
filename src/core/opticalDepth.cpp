#include "simprop/core/opticalDepth.h"

#include "simprop/crossSections/BreitWheeler.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace core {

double OpticalDepth::integrateOverPhField(double eGamma, double z, double mu) const {
  auto eps_th = 2. * pow2(SI::electronMassC2) / eGamma / (1. + z) / (1. - mu);
  auto epsMin = std::max(eps_th, m_ebl->getMinPhotonEnergy());
  auto epsMax = m_ebl->getMaxPhotonEnergy();
  return utils::QAGIntegration<double>(
      [this, eGamma, z, mu](double lnEps) {
        auto eps = std::exp(lnEps);
        auto n_gamma = m_ebl->density(eps, z);
        return eps * n_gamma * BreitWheeler::sigma(eGamma * (1. + z), eps, mu);
      },
      std::log(epsMin), std::log(epsMax), 1000, 1e-3);
}

double OpticalDepth::integrateOverAngle(double eGamma, double z) const {
  auto integrand = [this, eGamma, z](double mu) {
    return (1. - mu) * integrateOverPhField(eGamma, z, mu);
  };
  return utils::QAGIntegration<double>(integrand, -1, 1, 1000, 1e-3);
}

double OpticalDepth::get(double eGamma, double zSource) const {
  auto integrand = [this, eGamma](double z) {
    return m_cosmology->dtdz(z) * pow3(1. + z) * integrateOverAngle(eGamma, z);
  };
  auto value = utils::QAGIntegration<double>(integrand, 0., zSource, 1000, 1e-3);
  return 0.5 * SI::cLight * value;
}

}  // namespace core
}  // namespace simprop
