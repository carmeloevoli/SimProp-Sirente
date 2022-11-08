#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

#include "simprop/crossSections/PhotoPionProductionXsec.h"
#include "simprop/interactions/PhotoPionProduction.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/timer.h"

namespace simprop {
namespace losses {

PhotoPionContinuousLosses::PhotoPionContinuousLosses() : ContinuousLosses() {
  LOGD << "calling " << __func__ << " constructor";
  m_cmb = std::make_shared<photonfields::CMB>();
  m_sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  // const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  //  utils::Timer timer("caching table took : ");
  //  m_totalLosses.cacheTable(
  //      [&](double logGamma) {
  //        return m_inelasticity * pppcmb->rate(proton, std::pow(10., logGamma));
  //      },
  //      {9.95, 16.});
}

// double PhotoPionContinuousLosses::getInterpolated(double Gamma) const {
//   double b_l = 0;
//   const auto logGamma = std::log10(Gamma);
//   if (m_totalLosses.xIsInside(logGamma)) {
//     b_l = m_totalLosses.get(logGamma);
//   }
//   return b_l;
// }

double PhotoPionContinuousLosses::computeBetaComoving(double Gamma, double z) const {
  auto value = SI::cLight / 2. / pow2(Gamma);

  auto threshold = m_sigma->getPhotonEnergyThreshold();
  auto lnEpsPrimeMin = std::log(std::max(threshold, 2. * Gamma * m_cmb->getMinPhotonEnergy()));
  auto lnEpsPrimeMax = std::log(2. * Gamma * m_cmb->getMaxPhotonEnergy());

  value *= utils::simpsonIntegration<double>(
      [&](double logEpsPrime) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * m_sigma->getAtEpsPrime(epsPrime) *
               m_cmb->I_gamma(epsPrime / 2. / Gamma, z);
      },
      lnEpsPrimeMin, lnEpsPrimeMax, 300);

  return std::max(value, 0.);
}

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  auto A = (double)getPidNucleusMassNumber(pid);
  const double inelasticity = 0.15;
  return pow3(1. + z) * A * inelasticity * computeBetaComoving(Gamma * (1. + z), z);
}

}  // namespace losses
}  // namespace simprop
