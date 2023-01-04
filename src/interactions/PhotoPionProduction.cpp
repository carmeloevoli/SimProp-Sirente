#include "simprop/interactions/PhotoPionProduction.h"

#include <cmath>
#include <iostream>

#include "simprop/core/common.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

PhotoPionProduction::PhotoPionProduction(const std::shared_ptr<photonfields::PhotonField>& phField)
    : Interaction(phField) {
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoPionProduction::computeRateComoving(double Gamma, double z,
                                                const xsecs::CrossSection& xs) const {
  auto value = double(0);

  auto threshold = xs.getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_phField->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_phField->getMaxPhotonEnergy());
  if (lnEpsPrimeMax > lnEpsPrimeMin) {
    value = utils::simpsonIntegration<double>(
        [&](double lnEpsPrime) {
          auto epsPrime = std::exp(lnEpsPrime);
          auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
          return epsPrime * epsPrime * xs.getAtS(s) * m_phField->I_gamma(epsPrime / 2. / Gamma, z);
        },
        lnEpsPrimeMin, lnEpsPrimeMax, 300);
    value *= SI::cLight / 2. / pow2(Gamma);
  }
  return std::max(value, 0.);
}

double PhotoPionProduction::rate(PID pid, double Gamma, double z) const {
  const auto Z = getPidNucleusCharge(pid);
  const auto A = getPidNucleusMassNumber(pid);
  auto value = (double)Z * computeRateComoving(Gamma * (1. + z), z, m_xs_proton);
  if (A > Z) value += (double)(A - Z) * computeRateComoving(Gamma * (1. + z), z, m_xs_neutron);
  return pow3(1. + z) * value;
}

double sampleS(double r, double sMax, const xsecs::CrossSection& xs) {
  constexpr auto sThr = pow2(SI::protonMassC2 + SI::pionMassC2);
  if (sMax <= sThr) return 0;
  auto rPhiMax = r * xs.getPhiAtS(sMax);
  return utils::rootFinder<double>([&](double s) { return xs.getPhiAtS(s) - rPhiMax; }, sThr, sMax,
                                   1000, 1e-4);
}

double pickMinPhotonEnergy(double minPhotonEnergy, double nucleonEnergy) {
  const auto thresholdPhotonEnergy =
      (pow2(SI::pionMassC2) + 2. * SI::pionMassC2 * SI::protonMassC2) / 4. / nucleonEnergy;
  return std::max(minPhotonEnergy, thresholdPhotonEnergy);
}

double epsPdfIntegral(double photonEnergy, double nucleonEnergy, double z,
                      const xsecs::CrossSection& xs,
                      const std::shared_ptr<photonfields::PhotonField>& phField) {
  auto integrand = [&](double eps) {
    const auto s_max = pow2(SI::protonMassC2) + 4. * nucleonEnergy * eps;
    return phField->density(eps / (1. + z), z) / pow2(eps) * xs.getPhiAtS(s_max);
  };
  auto minPhEnergy = pickMinPhotonEnergy(phField->getMinPhotonEnergy(), nucleonEnergy);
  auto value = utils::QAGIntegration<double>(integrand, minPhEnergy, photonEnergy, 1000, 4e-4);
  // TODO try to improve the precision
  return value;
}

double sampleEps(double r, double nucleonEnergy, double z, const xsecs::CrossSection& xs,
                 const std::shared_ptr<photonfields::PhotonField>& phField) {
  auto minPhEnergy = pickMinPhotonEnergy(phField->getMinPhotonEnergy(), nucleonEnergy);
  auto maxPhotonEnergy = phField->getMaxPhotonEnergy();
  auto rIntegralMax = r * epsPdfIntegral(maxPhotonEnergy, nucleonEnergy, z, xs, phField);
  auto value = utils::rootFinder<double>(
      [&](double eps) { return epsPdfIntegral(eps, nucleonEnergy, z, xs, phField) - rIntegralMax; },
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

// Particle PhotoPionProduction::generatePion(PID nucleon, double nucleonEnergy, double z) const {
//   const xsecs::CrossSection& xs = (nucleon == proton) ? m_xs_proton : m_xs_neutron;
//   const auto photonEnergy = sampleEps(rng(), nucleonEnergy, z);
//   const auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
//   const auto s = sampleS(rng(), sMax);
//   const auto outPionEnergy = samplePionInelasticity(rng(), s) * nucleonEnergy;
//   const auto outPionCharge = samplePionCharge(rng(), (nucleon == neutron));
//   return {outPionCharge, zInteractionPoint, outPionEnergy / SI::pionMassC2, w};
//}

std::vector<Particle> PhotoPionProduction::finalState(const Particle& incomingParticle,
                                                      double zInteractionPoint,
                                                      RandomNumberGenerator& rng) const {
  const auto pid = incomingParticle.getPid();
  assert(pidIsNucleus(pid));
  const auto w = incomingParticle.getWeight();
  const auto Gamma = incomingParticle.getGamma();

  const auto nucleon = (pid == proton || pid == neutron) ? pid : pickNucleon(rng(), pid);
  const auto nucleonEnergy = Gamma * SI::protonMassC2;

  double s;
  if (nucleon == proton) {
    const auto photonEnergy =
        sampleEps(rng(), nucleonEnergy, zInteractionPoint, m_xs_proton, m_phField);
    const auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
    s = sampleS(rng(), sMax, m_xs_proton);
  } else {
    const auto photonEnergy =
        sampleEps(rng(), nucleonEnergy, zInteractionPoint, m_xs_neutron, m_phField);
    const auto sMax = pow2(SI::protonMassC2) + 4. * nucleonEnergy * photonEnergy;
    s = sampleS(rng(), sMax, m_xs_neutron);
  }

  const auto outPionEnergy = samplePionInelasticity(rng(), s) * nucleonEnergy;
  const auto outPionCharge = samplePionCharge(rng(), (nucleon == neutron));

  const auto outNucleonEnergy = nucleonEnergy - outPionEnergy;

  auto outPion = Particle(outPionCharge, zInteractionPoint, outPionEnergy / SI::pionMassC2, w);
  auto outNucleon = Particle(nucleon, zInteractionPoint, outNucleonEnergy / SI::protonMassC2, w);
  // auto outNuclues =
  //     Particle(removeNucleon(pid, nucleon), zInteractionPoint, Gamma, w);  // TODO finish this

  return {outNucleon, outPion};
}

}  // namespace interactions
}  // namespace simprop