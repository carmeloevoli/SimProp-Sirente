#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

#include "simprop/core/common.h"
#include "simprop/interactions/PhotoPionProduction.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/Dominguez2011PhotonField.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/timer.h"

namespace simprop {
namespace losses {

double averageAngle(double s) {
  if (s > 1e4 * SI::GeV2) return -1;
  const double b = 12. / SI::GeV2;
  double I_mu = utils::QAGIntegration<double>(
      [b, s](double mu) { return mu * std::exp(b * mu2t(mu, s)); }, -1., 1., 1000, 1e-5);
  double I = utils::QAGIntegration<double>([b, s](double mu) { return std::exp(b * mu2t(mu, s)); },
                                           -1., 1., 1000, 1e-5);
  return I_mu / I;
}

double inelasticity(double s) {
  if (s < 1.1519 * SI::GeV2) return 0.;
  auto sqrts = std::sqrt(s);
  auto E_pi = (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / sqrts;
  auto p_pi = std::sqrt((s - pow2(SI::protonMassC2 + SI::pionMassC2)) *
                        (s - pow2(SI::protonMassC2 - SI::pionMassC2))) /
              2. / sqrts;
  auto beta_pi = p_pi / E_pi;
  auto cosTheta_pi = averageAngle(s);
  return (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / s *
         (1. + beta_pi * cosTheta_pi);
}

double inelasticityFitted(double s) {
  if (s < 1.1519 * SI::GeV2) return 0.;
  // using https://en.wikipedia.org/wiki/Fano_resonance
  const double ds = s / SI::GeV2 - 1.15;
  const double q = 1.03;
  const double R = 0.80;
  auto value = pow2(0.5 * q * R + ds) / (pow2(0.5 * R) + pow2(ds));
  return 0.1014 * value / (s / SI::GeV2);
}

double inelasticityPoorApproximation(double s, double cosTheta_pi) {
  auto sqrts = std::sqrt(s);
  auto E_pi = (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / sqrts;
  auto p_pi = std::sqrt((s - pow2(SI::protonMassC2 + SI::pionMassC2)) *
                        (s - pow2(SI::protonMassC2 - SI::pionMassC2))) /
              2. / sqrts;
  auto beta_pi = p_pi / E_pi;
  return (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / s *
         (1. + beta_pi * cosTheta_pi);
}

PhotoPionContinuousLosses::PhotoPionContinuousLosses(
    const std::shared_ptr<photonfields::PhotonField>& photonField)
    : ContinuousLosses() {
  m_photonFields.push_back(photonField);
  LOGD << "calling " << __func__ << " constructor";
}

PhotoPionContinuousLosses::PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields)
    : ContinuousLosses(), m_photonFields(photonFields) {
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoPionContinuousLosses::computeBetaComoving(PID pid, double Gamma, double z) const {
  auto value = 0.;
  auto epsThr = m_xs.getPhotonEnergyThreshold();

  for (auto phField : m_photonFields) {
    auto lnEpsPrimeMin = std::log(std::max(epsThr, 2. * Gamma * phField->getMinPhotonEnergy()));
    auto lnEpsPrimeMax = std::log(2. * Gamma * phField->getMaxPhotonEnergy());
    if (lnEpsPrimeMax > lnEpsPrimeMin) {
      value += utils::simpsonIntegration<double>(
          [&](double lnEpsPrime) {
            auto epsPrime = std::exp(lnEpsPrime);
            auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
            return epsPrime * epsPrime * inelasticityFitted(s) * m_xs.getAtS(pid, s) *
                   phField->I_gamma(epsPrime / 2. / Gamma, z);
          },
          lnEpsPrimeMin, lnEpsPrimeMax, 500);
    }
  }
  return SI::cLight / 2. / pow2(Gamma) * std::max(value, 0.);
}

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  return pow3(1. + z) * computeBetaComoving(pid, Gamma * (1. + z), z);
}

}  // namespace losses
}  // namespace simprop
