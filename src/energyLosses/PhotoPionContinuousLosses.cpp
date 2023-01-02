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
  if (s < 1.16 * SI::GeV2) return 0.;
  if (s > 1e3 * SI::GeV2) return -1;
  const double b = 12. / SI::GeV2;
  double I_mu = utils::QAGIntegration<double>(
      [b, s](double mu) { return mu * std::exp(b * mu2t(mu, s)); }, -1., 1., 1000, 1e-8);
  double I = utils::QAGIntegration<double>([b, s](double mu) { return std::exp(b * mu2t(mu, s)); },
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
  m_xs_proton = std::make_shared<xsecs::PhotoPionProtonXsec>();
  m_xs_neutron = std::make_shared<xsecs::PhotoPionNeutronXsec>();
  LOGD << "calling " << __func__ << " constructor";
}

PhotoPionContinuousLosses::PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields)
    : ContinuousLosses(), m_photonFields(photonFields) {
  m_xs_proton = std::make_shared<xsecs::PhotoPionProtonXsec>();
  m_xs_neutron = std::make_shared<xsecs::PhotoPionNeutronXsec>();
  LOGD << "calling " << __func__ << " constructor";
}

double PhotoPionContinuousLosses::computeBetaComoving(
    double Gamma, const std::shared_ptr<xsecs::CrossSection>& xs) const {
  auto value = 0.;
  auto epsThr = xs->getPhotonEnergyThreshold();

  for (auto phField : m_photonFields) {
    auto lnEpsPrimeMin = std::log(std::max(epsThr, 2. * Gamma * phField->getMinPhotonEnergy()));
    auto lnEpsPrimeMax = std::log(2. * Gamma * phField->getMaxPhotonEnergy());
    if (lnEpsPrimeMax > lnEpsPrimeMin) {
      value += utils::simpsonIntegration<double>(
          [&](double lnEpsPrime) {
            auto epsPrime = std::exp(lnEpsPrime);
            auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
            return epsPrime * epsPrime * inelasticity(s) * xs->getAtS(s) *
                   phField->I_gamma(epsPrime / 2. / Gamma);
          },
          lnEpsPrimeMin, lnEpsPrimeMax, 500);
    }
  }
  return SI::cLight / 2. / pow2(Gamma) * std::max(value, 0.);
}

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  auto Z = getPidNucleusCharge(pid);
  auto A = getPidNucleusMassNumber(pid);
  auto value = (double)Z * computeBetaComoving(Gamma * (1. + z), m_xs_proton);
  if (A > Z) value += (double)(A - Z) * computeBetaComoving(Gamma * (1. + z), m_xs_neutron);
  return pow3(1. + z) * value;
}

}  // namespace losses
}  // namespace simprop
