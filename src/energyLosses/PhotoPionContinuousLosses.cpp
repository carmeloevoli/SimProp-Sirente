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
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  utils::Timer timer("caching table took : ");
  m_totalLosses.cacheTable(
      [&](double logGamma) {
        return m_inelasticity * pppcmb->rate(proton, std::pow(10., logGamma));
      },
      {9.95, 16.});
}

double PhotoPionContinuousLosses::getInterpolated(double Gamma) const {
  double b_l = 0;
  const auto logGamma = std::log10(Gamma);
  if (m_totalLosses.xIsInside(logGamma)) {
    b_l = m_totalLosses.get(logGamma);
  }
  return b_l;
}

double PhotoPionContinuousLosses::beta(PID pid, double Gamma, double z) const {
  double b_l = getInterpolated(Gamma * (1. + z));
  if (b_l > 0.) {
    b_l *= pow3(1. + z);
    const double Z = (double)getPidNucleusCharge(pid);
    const double A = (double)getPidNucleusMassNumber(pid);
    b_l *= pow2(Z) / A;
  }
  return std::max(b_l, 0.);
}

}  // namespace losses
}  // namespace simprop
