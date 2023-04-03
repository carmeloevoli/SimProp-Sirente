#include "simprop/interactions/PhotoPionProduction.h"

#include <cmath>
#include <iostream>

#include "simprop/core/common.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

double PhotoPionProduction::sampleS(double r, PID nucleon, double sMax) const {
  constexpr auto sThr = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (sMax <= sThr) return 0;
  auto rPhiMax = r * m_xs.getPhiAtS(nucleon, sMax);
  return utils::rootFinder<double>([&](double s) { return m_xs.getPhiAtS(nucleon, s) - rPhiMax; },
                                   sThr, sMax, 1000, 1e-4);
}

double pickMinPhotonEnergy(double minPhotonEnergy, double nucleonEnergy) {
  const auto thresholdPhotonEnergy =
      (pow2(SI::pionMassC2) + 2. * SI::pionMassC2 * SI::protonMassC2) / 4. / nucleonEnergy;
  return std::max(minPhotonEnergy, thresholdPhotonEnergy);
}

double PhotoPionProduction::epsPdfIntegral(double photonEnergy, PID nucleon, double nucleonEnergy,
                                           double z) const {
  auto integrand = [&](double eps) {
    const auto s_max = pow2(SI::protonMassC2) + 4. * nucleonEnergy * eps;
    return m_phField->density(eps / (1. + z), z) / pow2(eps) * m_xs.getPhiAtS(nucleon, s_max);
  };
  auto minPhEnergy = pickMinPhotonEnergy(m_phField->getMinPhotonEnergy(), nucleonEnergy);
  auto value = utils::QAGIntegration<double>(integrand, minPhEnergy, photonEnergy, 1000, 4e-4);
  // TODO how to improve the precision here?
  return value;
}

double PhotoPionProduction::sampleEps(double r, PID nucleon, double nucleonEnergy, double z) const {
  auto minPhEnergy = pickMinPhotonEnergy(m_phField->getMinPhotonEnergy(), nucleonEnergy);
  auto maxPhotonEnergy = m_phField->getMaxPhotonEnergy();
  auto rIntegralMax = r * epsPdfIntegral(maxPhotonEnergy, nucleon, nucleonEnergy, z);
  auto value = utils::rootFinder<double>(
      [&](double eps) { return epsPdfIntegral(eps, nucleon, nucleonEnergy, z) - rIntegralMax; },
      minPhEnergy, maxPhotonEnergy, 1000, 1e-3);
  return value;
}

double angleCoMintegral(double mu, double s) {
  constexpr auto b = 12. / SI::GeV2;
  auto integrand = [b, s](double mu) { return std::exp(b * mu2t(mu, s)); };
  auto value = utils::QAGIntegration<double>(integrand, -1., mu, 1000, 1e-5);
  return value;
}

double sampleAngleCoM(double r, double s) {
  auto rIntegralMax = r * angleCoMintegral(1., s);
  auto value = utils::rootFinder<double>(
      [s, rIntegralMax](double mu) { return angleCoMintegral(mu, s) - rIntegralMax; }, -1., 1.,
      1000, 1e-3);
  return value;
}

double samplePionInelasticity(double r, double s) {
  auto sqrt_s = std::sqrt(s);
  auto E_star = 0.5 * (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / sqrt_s;
  auto p_star = 0.5 / sqrt_s;
  p_star *= std::sqrt((s - pow2(SI::pionMassC2 + SI::protonMassC2)) *
                      (s - pow2(SI::pionMassC2 - SI::protonMassC2)));
  // auto mu_star = 2. * r - 1.;
  auto mu_star = sampleAngleCoM(r, s);
  return 1. / sqrt_s * (E_star + p_star * mu_star);
}

PID samplePionCharge(double r, bool isNeutron) {
  if (r < 2. / 3.)
    return pionNeutral;
  else
    return (isNeutron) ? pionMinus : pionPlus;
}

PID pickNucleon(double r, PID pid) {
  auto Z = (double)getPidNucleusCharge(pid);
  auto A = (double)getPidNucleusMassNumber(pid);
  if (r < Z / A)
    return proton;
  else
    return neutron;
}

PhotoPionProduction::PhotoPionProduction(const std::shared_ptr<photonfields::PhotonField>& phField)
    : Interaction(phField) {
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoPionProduction::rate(PID pid, double Gamma, double z) const {
  auto value = double(0);

  auto threshold = m_xs.getEpsPrimeThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_phField->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_phField->getMaxPhotonEnergy());
  if (lnEpsPrimeMax > lnEpsPrimeMin) {
    value = utils::simpsonIntegration<double>(
        [&](double lnEpsPrime) {
          auto epsPrime = std::exp(lnEpsPrime);
          return epsPrime * epsPrime * m_xs.getAtEpsPrime(pid, epsPrime) *
                 m_phField->I_gamma(epsPrime / 2. / Gamma, z);
        },
        lnEpsPrimeMin, lnEpsPrimeMax, 500);
    value *= SI::cLight / 2. / pow2(Gamma);
  }
  return std::max(value, 0.);
}

std::vector<Particle> PhotoPionProduction::finalState(const Particle& incomingParticle,
                                                      double zInteractionPoint,
                                                      RandomNumberGenerator& rng) const {
  const auto pid = incomingParticle.getPid();
  assert(pidIsNucleus(pid));
  const auto w = incomingParticle.getWeight();
  const auto Gamma = incomingParticle.getGamma();

  const auto nucleon = pidIsNucleon(pid) ? pid : pickNucleon(rng(), pid);
  const auto nucleonEnergy = Gamma * SI::protonMassC2;

  const auto photonEnergy = sampleEps(rng(), nucleon, nucleonEnergy, zInteractionPoint);
  const auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
  const auto s = sampleS(rng(), nucleon, sMax);

  const auto outPionEnergy = samplePionInelasticity(rng(), s) * nucleonEnergy;
  const auto outPionCharge = samplePionCharge(rng(), (nucleon == neutron));

  const auto outNucleonEnergy = nucleonEnergy - outPionEnergy;

  auto outPion = Particle(outPionCharge, zInteractionPoint, outPionEnergy / SI::pionMassC2, w);
  auto outNucleon = Particle(proton, zInteractionPoint, outNucleonEnergy / SI::protonMassC2, w);

  if (pidIsNucleon(pid)) {
    return {outNucleon, outPion};
  } else {
    auto outNuclues = Particle(removeNucleon(pid, nucleon), zInteractionPoint, Gamma, w);
    return {outNuclues, outNucleon, outPion};
  }
}

}  // namespace interactions
}  // namespace simprop