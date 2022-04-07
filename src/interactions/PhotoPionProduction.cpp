#include "simprop/interactions/PhotoPionProduction.h"

#include <cmath>
#include <iostream>

#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double PhotoPionProduction::computeRateComoving(PID pid, double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  auto threshold = m_sigma->getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_ebl->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_ebl->getMaxPhotonEnergy());

  value *= utils::simpsonIntegration<double>(
      [&](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * m_sigma->getAtEpsPrime(proton, epsPrime) *
               m_ebl->I_gamma(epsPrime / 2. / Gamma, z);
      },
      lnEpsPrimeMin, lnEpsPrimeMax, 300);

  return std::max(value, 0.);
}

double PhotoPionProduction::rate(PID pid, double Gamma, double z) const {
  return pow3(1. + z) * computeRateComoving(pid, Gamma * (1. + z), z);
}

double PhotoPionProduction::rateCrpropaComparison(PID pid, double Gamma, double z) const {
  return computeRateComoving(pid, Gamma, z);
}

double PhotoPionProduction::phi(double s) const {
  auto integrand = [&](double s) {
    return (s - pow2(SI::protonMassC2)) * m_sigma->getAtS(proton, s);
  };
  auto value = utils::QAGIntegration<double>(integrand, m_sThreshold, s, 1000, 0.5e-3);
  return value;
}

double PhotoPionProduction::sample_s(double r, double sMax) const {
  if (sMax <= m_sThreshold) return 0;
  auto rPhiMax = r * phi(sMax);
  return utils::rootFinder<double>([&](double s) { return phi(s) - rPhiMax; }, m_sThreshold, sMax,
                                   1000, 1e-4);
}

double PhotoPionProduction::epsPdfIntegral(double photonEnergy, double E, double z) const {
  auto integrand = [&](double eps) {
    auto s_max = pow2(SI::protonMassC2) + 4. * E * eps;
    return m_ebl->density(eps) / pow2(eps) * phi(s_max);
  };
  auto value = utils::QAGIntegration<double>(integrand, m_ebl->getMinPhotonEnergy(), photonEnergy,
                                             1000, 0.5e-3);
  return value;
}

double PhotoPionProduction::sample_eps(double r, double E, double z) const {
  auto rIntegralMax = r * epsPdfIntegral(m_ebl->getMaxPhotonEnergy(), E, z);
  auto value = utils::rootFinder<double>(
      [&](double eps) { return epsPdfIntegral(eps, E, z) - rIntegralMax; },
      m_ebl->getMinPhotonEnergy(), m_ebl->getMaxPhotonEnergy(), 1000, 1e-4);
  return value;
}

double PhotoPionProduction::samplePionEnergy(double E, double z, RandomNumberGenerator& rng) const {
  auto photonEnergy = sample_eps(rng(), E, z);
  auto sMax = pow2(SI::protonMassC2) + 4. * E * photonEnergy;
  auto s = sample_s(rng(), sMax);
  auto sqrt_s = std::sqrt(s);
  auto E_star = 0.5 * (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / sqrt_s;
  auto p_star = 0.5 *
                std::sqrt((s - pow2(SI::pionMassC2 + SI::protonMassC2)) *
                          (s - pow2(SI::pionMassC2 - SI::protonMassC2))) /
                sqrt_s;
  auto mu_star = rng.uniform(-1, 1);
  return E / sqrt_s * (E_star + p_star * mu_star);
}

}  // namespace interactions
}  // namespace simprop