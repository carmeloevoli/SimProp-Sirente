#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

#include "simprop/crossSections/PhotoPionProductionXsec.h"
#include "simprop/interactions/PhotoPionProduction.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/Dominguez2011PhotonField.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/timer.h"

namespace simprop {
namespace losses {

double t(double mu, double s) {
  const auto mp2 = pow2(SI::protonMassC2);
  const auto mpi2 = pow2(SI::pionMassC2);
  const auto sqrts = std::sqrt(s);
  double value = mp2;
  value -= (s + mp2) * (s + mp2 - mpi2) / 2. / s;
  value -= (s - mp2) / sqrts * std::sqrt(pow2((s + mp2 - mpi2) / 2. / sqrts) - mp2) * mu;
  return value;
}

double averageAngle(double s) {
  if (s < 1.16 * SI::GeV2) return 0.;
  if (s > 1e3 * SI::GeV2) return -1;
  const double b = 12. / SI::GeV2;
  double I_mu = utils::QAGIntegration<double>(
      [b, s](double mu) { return mu * std::exp(b * t(mu, s)); }, -1., 1., 1000, 1e-8);
  double I = utils::QAGIntegration<double>([b, s](double mu) { return std::exp(b * t(mu, s)); },
                                           -1., 1., 1000, 1e-8);
  return I_mu / I;
}

double inelasticity(double s) {
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

double inelasticityPoorApproximation(double s) {
  auto sqrts = std::sqrt(s);
  auto E_pi = (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / sqrts;
  auto p_pi = std::sqrt((s - pow2(SI::protonMassC2 + SI::pionMassC2)) *
                        (s - pow2(SI::protonMassC2 - SI::pionMassC2))) /
              2. / sqrts;
  auto beta_pi = p_pi / E_pi;
  auto cosTheta_pi = -1.;
  return (s - pow2(SI::protonMassC2) + pow2(SI::pionMassC2)) / 2. / s *
         (1. + beta_pi * cosTheta_pi);
}

PhotoPionContinuousLosses::PhotoPionContinuousLosses(
    const std::shared_ptr<photonfields::PhotonField>& photonField)
    : ContinuousLosses() {
  m_photonFields.push_back(photonField);
  m_sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();

  LOGD << "calling " << __func__ << " constructor";
}

PhotoPionContinuousLosses::PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields)
    : ContinuousLosses(), m_photonFields(photonFields) {
  m_sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  LOGD << "calling " << __func__ << " constructor";
}

// const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
//  utils::Timer timer("caching table took : ");
//  m_totalLosses.cacheTable(
//      [&](double logGamma) {
//        return m_inelasticity * pppcmb->rate(proton, std::pow(10., logGamma));
//      },
//      {9.95, 16.});

// double PhotoPionContinuousLosses::getInterpolated(double Gamma) const {
//   double b_l = 0;
//   const auto logGamma = std::log10(Gamma);
//   if (m_totalLosses.xIsInside(logGamma)) {
//     b_l = m_totalLosses.get(logGamma);
//   }
//   return b_l;
// }

double PhotoPionContinuousLosses::computeBetaComoving(double Gamma) const {
  auto value = 0.;
  auto threshold = m_sigma->getPhotonEnergyThreshold();
  for (auto phField : m_photonFields) {
    auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * phField->getMinPhotonEnergy()));
    auto lnEpsPrimeMax = std::log(2. * Gamma * phField->getMaxPhotonEnergy());

    value += utils::simpsonIntegration<double>(
        [&](double logEpsPrime) {
          auto epsPrime = std::exp(logEpsPrime);
          auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
          return epsPrime * epsPrime * inelasticity(s) * m_sigma->getAtEpsPrime(epsPrime) *
                 phField->I_gamma(epsPrime / 2. / Gamma);
        },
        lnEpsPrimeMin, lnEpsPrimeMax, 500);
  }
  return SI::cLight / 2. / pow2(Gamma) * std::max(value, 0.);
}

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  auto A = (double)getPidNucleusMassNumber(pid);
  return pow3(1. + z) * A * computeBetaComoving(Gamma * (1. + z));
}

}  // namespace losses
}  // namespace simprop
