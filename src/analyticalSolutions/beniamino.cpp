#include "simprop/analyticalSolutions/beniamino.h"

#include <limits>
#include <numeric>

#include "simprop/utils/logging.h"

namespace simprop {
namespace solutions {

#define INTSTEPS 1000
#define VERYLARGEENERGY (1e25 * SI::eV)
#define VERYLARGEJACOBIAN (1e6)

Beniamino::Beniamino(const SourceParams& params, const std::shared_ptr<cosmo::Cosmology>& cosmology,
                     const std::vector<std::shared_ptr<losses::ContinuousLosses>>& losses)
    : m_cosmology(cosmology), m_losses(losses) {
  m_injSlope = params.injSlope;
  m_evolutionIndex = params.evolutionIndex;
  m_expCutoff = params.expCutoff;
  LOGD << "calling " << __func__ << " constructor";
}

Beniamino& Beniamino::doCaching() {
  m_lossesLookup.cacheTable(
      [this](double lnE) {
        const auto Gamma = std::exp(lnE) / SI::protonMassC2;
        return std::accumulate(
            m_losses.begin(), m_losses.end(), 0.,
            [Gamma](double sum, const std::shared_ptr<losses::ContinuousLosses>& l) {
              return sum + l->beta(proton, Gamma);
            });
      },
      {std::log(1e16 * SI::eV), std::log(VERYLARGEENERGY)});
  m_doCaching = true;
  return *this;
}

double Beniamino::beta(double E) const {
  if (E < 1e16 * SI::eV || E > VERYLARGEENERGY) return 0;
  if (m_doCaching)
    return m_lossesLookup.get(std::log(E));
  else {
    const auto Gamma = E / SI::protonMassC2;
    return std::accumulate(m_losses.begin(), m_losses.end(), 0.,
                           [Gamma](double sum, const std::shared_ptr<losses::ContinuousLosses>& l) {
                             return sum + l->beta(proton, Gamma);
                           });
  }
}

double Beniamino::generationEnergy(double E, double zNow, double zMax, double relError) const {
  assert(zMax >= zNow);
  assert(E > 0);
  auto dEdz = [this](double z, double E_g) {
    auto dtdz = m_cosmology->dtdz(z);
    auto E = std::min(E_g * (1. + z), VERYLARGEENERGY);
    return E_g * (1. / (1. + z) + dtdz * pow3(1. + z) * beta(E));
  };
  auto value = utils::odeiv<double>(dEdz, E, zNow, zMax, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEENERGY);
}

double Beniamino::dbdE(double E) const {
  if (E < 1e16 * SI::eV || E > VERYLARGEENERGY) return 0;
  auto factor = utils::deriv<double>(
      [this](double lnx) {
        auto x = std::exp(lnx);
        auto betax = std::max(beta(x), 1e-50 / SI::year);
        return std::log(x * betax);
      },
      std::log(E), 1e-2);
  return beta(E) * factor;
}

double Beniamino::dilationFactor(double E, double zNow, double zMax, double relError) const {
  assert(zMax >= zNow);
  assert(E > 0);
  if (generationEnergy(E, zNow, zMax, 1e-2) > 0.2 * VERYLARGEENERGY) return VERYLARGEJACOBIAN;
  auto dydz = [this, E, zNow](double z, double y) {
    auto dtdz = m_cosmology->dtdz(z);
    auto E_g = generationEnergy(E, zNow, z, 1e-5);
    auto E_prime = std::min(E_g * (1. + z), VERYLARGEENERGY);
    auto dbdE = std::max(this->dbdE(E_prime), 0.);
    return y * (1. / (1. + z) + dtdz * pow3(1. + z) * dbdE);
  };
  auto value = utils::odeiv<double>(dydz, 1., zNow, zMax, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEJACOBIAN);
}

double Beniamino::computeFlux(double E, double zObs, double zMax, double relError) const {
  const auto K = (m_injSlope - 2.) / pow2(m_minEnergy);
  const auto L_0 = m_sourceEmissivity;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0;
  auto integrand = [this, E, zObs](double z) {
    const auto E_g = generationEnergy(E, zObs, z, 1e-4);
    if (E_g > m_maxEnergy) return 0.;
    auto dEgdE = dilationFactor(E, zObs, z, 1e-4);
    auto inj = std::pow(E_g / m_minEnergy, -m_injSlope);
    if (m_expCutoff > 0.) inj *= std::exp(-E_g / m_expCutoff);
    auto sourceEvolution = std::pow(1. + z, m_evolutionIndex);
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz * sourceEvolution * inj * dEgdE;
  };
  auto I = utils::QAGIntegration<double>(integrand, zObs, zMax, INTSTEPS, relError);

  return factor * I;
}

// double Beniamino::computeFluxUnm(double E, double zMax, double relError) const {
//   const auto K = m_slope - 2;
//   const auto L_0 = m_sourceEmissivity;
//   const auto factor = SI::cLight / 4. / M_PI * K * L_0 * std::pow(E / E_0, -m_slope);
//   auto integrand = [&](double z) { return m_cosmology->dtdz(z) * std::pow(1. + z, -m_slope
//   + 1.);
//   }; auto I = utils::QAGIntegration<double>(integrand, 0., zMax, INTSTEPS, relError);

//   return factor * I;
// }

}  // namespace solutions
}  // namespace simprop