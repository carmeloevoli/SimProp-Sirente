#include "simprop/energyLosses/PairProductionLosses.h"

#include <algorithm>
#include <array>

#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

double phi(double k) {
  auto _c = std::array<double, 4>{0.8048, 0.1459, 1.137e-3, -3.879e-6};
  auto _d = std::array<double, 4>{-86.07, 50.96, -14.45, 8 / 3.};
  auto _f = std::array<double, 3>{2.910, 78.35, 1837};
  if (k < 25) {
    double value = 0;
    for (size_t i = 1; i < 5; ++i) value += _c.at(i - 1) * std::pow(k - 2, (double)i);
    return M_PI / 12. * pow4(k - 2) / (1. + value);
  } else {
    // double phi_ = 0;
    // for (size_t i = 0; i < 4; ++i) phi_ += _d.at(i) * std::log(k, (double)i);
    // phi_ *= k;
    // double value = 0;
    // for (size_t i = 1; i < 4; ++i) value += _f.at(i) * std::pow(k, -(double)i);
    return 0;
  }
}

double PairProductionLosses::dlnE_dt(PID pid, double E, double z) const {
  auto Gamma = getGamma(pid, E);
  const auto cmb = photonfields::CMB();
  const auto epsmin = cmb.getMinPhotonEnergy();
  const auto epsmax = cmb.getMaxPhotonEnergy();

  const auto lkmin = 0;  // std::log(std::max(2, 2 * Gamma * epsmin));
  const auto lkmax = std::log(2 * Gamma * epsmax);

  auto I = utils::simpsonIntegration<double>(
      [Gamma, &cmb](double logk) {
        auto epsPrime = std::exp(logEpsPrime);
        return epsPrime * epsPrime * sigma.get(proton, epsPrime) *
               cmb.density(epsPrime / 2. / Gamma);
      },
      lkmin, lkmax, 100);

  const auto factor = SI::alpha * pow2(SI::electronRadius) * (SI::electronMass / SI::protonMass);
  return factor;
}

//   const auto redshiftedEnergy = E * (1. + z);
//   double b_l = getInterpolated(redshiftedEnergy);
//   if (b_l > 0.) {
//     b_l *= pow3(1. + z);
//     const double Z = (double)getNucleusChargeNumber(pid);
//     const double A = (double)getNucleusChargeNumber(pid);
//     b_l *= pow2(Z) / A;
//   }
// return 0;

double PairProductionLosses::dlnE_dz(PID pid, double E, double z) const {
  auto b_l = dlnE_dt(pid, E, z);
  return (b_l > 0.) ? b_l * m_cosmology.dtdz(z) : 0.;
}

}  // namespace losses
}  // namespace simprop
