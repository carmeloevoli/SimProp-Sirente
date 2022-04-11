#include "simprop/interactions/PhotoPionProduction.h"

#include <cmath>
#include <iostream>

#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double PhotoPionProduction::computeRateComoving(double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  auto threshold = m_sigma->getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_phField->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_phField->getMaxPhotonEnergy());

  value *= utils::simpsonIntegration<double>(
      [&](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * m_sigma->getAtEpsPrime(proton, epsPrime) *
               m_phField->I_gamma(epsPrime / 2. / Gamma, z);
      },
      lnEpsPrimeMin, lnEpsPrimeMax, 300);

  return std::max(value, 0.);
}

double PhotoPionProduction::rate(PID pid, double Gamma, double z) const {
  auto A = (double)getNucleusMassNumber(pid);
  return pow3(1. + z) * A * computeRateComoving(Gamma * (1. + z), z);
}

double PhotoPionProduction::phi(double s) const {
  auto integrand = [&](double s) {
    return (s - pow2(SI::protonMassC2)) * m_sigma->getAtS(proton, s);
  };
  auto value = utils::QAGIntegration<double>(integrand, m_sThreshold, s, 1000, 0.5e-3);
  return value;
}

double PhotoPionProduction::sampleS(RndUnifNumber r, double sMax) const {
  if (sMax <= m_sThreshold) return 0;
  auto rPhiMax = r.get() * phi(sMax);
  return utils::rootFinder<double>([&](double s) { return phi(s) - rPhiMax; }, m_sThreshold, sMax,
                                   1000, 1e-4);
}

double pickMinPhotonEnergy(double minPhotonEnergy, double nucleonEnergy) {
  const auto thresholdPhotonEnergy =
      (pow2(SI::pionMassC2) + 2. * SI::pionMassC2 * SI::protonMassC2) / 4. / nucleonEnergy;
  return std::max(minPhotonEnergy, thresholdPhotonEnergy);
}

double PhotoPionProduction::epsPdfIntegral(double photonEnergy, double nucleonEnergy,
                                           double z) const {
  auto integrand = [&](double eps) {
    const auto s_max = pow2(SI::protonMassC2) + 4. * nucleonEnergy * eps;
    return m_phField->density(eps) / pow2(eps) * phi(s_max);  // TODO evolution in z?
  };
  auto minPhEnergy = pickMinPhotonEnergy(m_phField->getMinPhotonEnergy(), nucleonEnergy);
  auto value = utils::QAGIntegration<double>(integrand, minPhEnergy, photonEnergy, 1000,
                                             1e-3);  // TODO why improve this?
  return value;
}

double PhotoPionProduction::sampleEps(RndUnifNumber r, double nucleonEnergy, double z) const {
  auto minPhEnergy = pickMinPhotonEnergy(m_phField->getMinPhotonEnergy(), nucleonEnergy);
  auto maxPhotonEnergy = m_phField->getMaxPhotonEnergy();
  auto rIntegralMax = r.get() * epsPdfIntegral(maxPhotonEnergy, nucleonEnergy, z);
  auto value = utils::rootFinder<double>(
      [&](double eps) { return epsPdfIntegral(eps, nucleonEnergy, z) - rIntegralMax; }, minPhEnergy,
      maxPhotonEnergy, 1000, 1e-3);
  return value;
}

double PhotoPionProduction::samplePionInelasticity(RndUnifNumber r, double s) const {
  // auto photonEnergy = sample_eps(rng(), nucleonEnergy, z);
  // auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
  // auto s = sample_s(rng(), sMax);
  auto sqrt_s = std::sqrt(s);
  auto E_star = 0.5 * (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / sqrt_s;
  auto p_star = 0.5 *
                std::sqrt((s - pow2(SI::pionMassC2 + SI::protonMassC2)) *
                          (s - pow2(SI::pionMassC2 - SI::protonMassC2))) /
                sqrt_s;
  auto mu_star = 2. * r.get() - 1.;
  return 1. / sqrt_s * (E_star + p_star * mu_star);
}

PID PhotoPionProduction::samplePionCharge(RndUnifNumber r, bool isNeutron) const {
  if (r.get() < 2. / 3.)
    return pionNeutral;
  else
    return (isNeutron) ? pionMinus : pionPlus;
}

std::vector<Particle> PhotoPionProduction::finalState(const Particle& incomingParticle,
                                                      double zInteractionPoint,
                                                      RandomNumberGenerator& rng) const {
  const auto pid = incomingParticle.getPid();
  if (pid == proton || pid == neutron) {
    const auto nucleonEnergy = incomingParticle.getGamma() * SI::protonMassC2;
    const auto photonEnergy =
        sampleEps(rng(), nucleonEnergy, zInteractionPoint);  // TODO implement z-evolution
    const auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
    const auto s = sampleS((rng()), sMax);
    const auto outPionEnergy = samplePionInelasticity(rng(), s) * nucleonEnergy;
    const auto outPionCharge = samplePionCharge(rng(), (pid == neutron));
    const auto outNucleonEnergy = nucleonEnergy - outPionEnergy;
    auto outPion = Particle(outPionCharge, zInteractionPoint, outPionEnergy / SI::pionMassC2);
    auto outNucleon = Particle(proton, zInteractionPoint, outNucleonEnergy / SI::protonMassC2);
    return {outNucleon, outPion};
  }
  auto p = Particle{incomingParticle};
  return {p};
}

}  // namespace interactions
}  // namespace simprop