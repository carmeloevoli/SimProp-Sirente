// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

#include "simprop/core/common.h"
#include "simprop/interactions/PhotoPionProduction.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/Dominguez2011PhotonField.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/timer.h"

namespace simprop {
namespace losses {

double inelasticity(double epsPrime) {
  static const double Y0 = 0.47;
  static const double b = 6e9 * SI::eV;
  static const double d = 0.33;
  static const double s = 0.15;
  const auto x = epsPrime / b;
  return Y0 * std::pow(x, d) / std::pow(1. + std::pow(x, d / s), s);
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

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  auto value = 0.;
  auto epsThr = m_xs.getEpsPrimeThreshold();

  for (auto phField : m_photonFields) {
    auto lnEpsPrimeMin = std::log(std::max(epsThr, 2. * Gamma * phField->getMinPhotonEnergy()));
    auto lnEpsPrimeMax = std::log(2. * Gamma * phField->getMaxPhotonEnergy());
    if (lnEpsPrimeMax > lnEpsPrimeMin) {
      value += utils::RombergIntegration<double>(
          [&](double lnEpsPrime) {
            auto epsPrime = std::exp(lnEpsPrime);
            return epsPrime * epsPrime * m_xs.getAtEpsPrime(pid, epsPrime) *
                   inelasticity(epsPrime) * phField->I_gamma(epsPrime / 2. / Gamma, z);
          },
          lnEpsPrimeMin, lnEpsPrimeMax, 9, 1e-3);
    }
  }
  return SI::cLight / 2. / pow2(Gamma) * std::max(value, 0.);
}

}  // namespace losses
}  // namespace simprop
