#include "simprop/energyLosses/PairProductionLosses.h"

#include <array>
#include <cmath>

#include "simprop/core/units.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

double sum_c(double k) {
  static auto _c = std::array<double, 4>{0.8048, 0.1459, 1.137e-3, -3.879e-6};
  double value = 0;
  for (size_t i = 1; i <= 4; ++i) value += _c[i - 1] * std::pow(k - 2, (double)i);
  return value;
}

double sum_d(double k) {
  static auto _d = std::array<double, 4>{-86.07, 50.96, -14.45, 8 / 3.};
  auto lnk = std::log(k);
  double value = 0;
  for (size_t i = 0; i <= 3; ++i) value += _d[i] * std::pow(lnk, (double)i);
  return value;
}

double sum_f(double k) {
  static auto _f = std::array<double, 3>{2.910, 78.35, 1837};
  double value = 0;
  for (size_t i = 1; i <= 3; ++i) value += _f[i - 1] / std::pow(k, (double)i);
  return value;
}

double phi(double k) {
  if (k < 2.)
    return 0;
  else if (k < 25.) {
    return M_PI / 12. * pow4(k - 2) / (1. + sum_c(k));
  } else {
    return (k * sum_d(k)) / (1. - sum_f(k));
  }
}

PairProductionLosses::PairProductionLosses(
    const std::shared_ptr<photonfields::PhotonField>& photonField)
    : ContinuousLosses() {
  m_photonFields.push_back(photonField);
  LOGD << "calling " << __func__ << " constructor";
}

PairProductionLosses::PairProductionLosses(const photonfields::PhotonFields& photonFields)
    : ContinuousLosses(), m_photonFields(photonFields) {
  LOGD << "calling " << __func__ << " constructor";
}

PairProductionLosses& PairProductionLosses::doCaching() {
  m_betaProtons.cacheTable(
      [this](double lnGamma, double z) {
        auto Gamma = std::exp(lnGamma);
        return computeProtonBeta(Gamma, z);
      },
      {std::log(1e7), std::log(1e14)}, {0., 10.});
  m_doCaching = true;
  return *this;
}

double PairProductionLosses::computeProtonBeta(double Gamma, double z) const {
  auto TwoGamma_mec2 = 2. * Gamma / SI::electronMassC2;
  double value = 0;
  for (auto phField : m_photonFields) {
    const auto epsmin = phField->getMinPhotonEnergy();
    const auto epsmax = phField->getMaxPhotonEnergy();
    const auto lkmin = std::log(TwoGamma_mec2 * epsmin);
    const auto lkmax = std::log(TwoGamma_mec2 * epsmax);
    value += utils::RombergIntegration<double>(
        [TwoGamma_mec2, phField, z](double lnk) {
          auto k = std::exp(lnk);
          return phi(k) / k * phField->density(k / TwoGamma_mec2, z);
        },
        lkmin, lkmax, 8, 1e-3);
  }
  constexpr auto factor = SI::alpha * pow2(SI::electronRadius) * SI::cLight * SI::electronMassC2 *
                          (SI::electronMass / SI::protonMass);
  return factor * value / Gamma;
}

double PairProductionLosses::beta(PID pid, double Gamma, double z) const {
  auto b_l = (m_doCaching) ? m_betaProtons.get(std::log(Gamma), z) : computeProtonBeta(Gamma, z);
  auto Z = (double)getPidNucleusCharge(pid);
  auto A = (double)getPidNucleusMassNumber(pid);
  b_l *= pow2(Z) / A;
  return std::max(b_l, 0.);
}

}  // namespace losses
}  // namespace simprop
