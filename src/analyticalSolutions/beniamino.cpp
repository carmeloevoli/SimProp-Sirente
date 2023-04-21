#include "simprop/analyticalSolutions/beniamino.h"

#include <limits>

#include "simprop/utils/logging.h"

namespace simprop {
namespace solutions {

#define INTSTEPS 1000

#define VERYLARGEENERGY (1e25 * SI::eV)
#define VERYLARGEJACOBIAN (1e6)

Beniamino::Beniamino(bool doPhotoPion) {
  m_cosmology = std::make_shared<cosmo::Cosmology>();
  {
    auto cmb = std::make_shared<photonfields::CMB>();
    auto pair = std::make_shared<losses::PairProductionLosses>(cmb);
    auto pion = std::make_shared<losses::PhotoPionContinuousLosses>(cmb);
    m_losses.cacheTable(
        [doPhotoPion, pair, pion](double lnE) {
          auto Gamma = std::exp(lnE) / SI::protonMassC2;
          return pair->beta(proton, Gamma) + ((doPhotoPion) ? pion->beta(proton, Gamma) : 0.);
        },
        {std::log(1e16 * SI::eV), std::log(VERYLARGEENERGY)});
  }
  LOGD << "calling " << __func__ << " constructor";
}

Beniamino::Beniamino(bool doPhotoPion, BeniaminoParams params) : Beniamino(doPhotoPion) {
  m_slope = params.slope;
  m_sourceEvolution = params.sourceEvolution;
  m_sourceCutoff = params.sourceCutoff;
  LOGD << "calling " << __func__ << " constructor";
}

double Beniamino::dbdE(double E) const {
  if (E < 1e17 * SI::eV || E > VERYLARGEENERGY) return 0;
  auto dbetadE =
      utils::deriv5pt<double>([this](double x) { return m_losses.get(std::log(x)); }, E, 1e-1 * E);

  return m_losses.get(std::log(E)) + E * dbetadE;
}

double Beniamino::beta(double E) const { return m_losses.get(std::log(E)); }

double Beniamino::generationEnergy(double E, double zNow, double zMax, double relError) const {
  assert(zMax >= zNow);
  assert(E > 0);
  auto dEdz = [this, zNow](double z, double E_g) {
    auto E = std::min(E_g * (1. + z) / (1. + zNow), VERYLARGEENERGY);
    auto beta = m_losses.get(std::log(E));
    auto dtdz = m_cosmology->dtdz(z);
    return E_g * (1. / (1. + z) + dtdz * pow3(1. + z) * beta);
  };
  auto value = utils::odeiv<double>(dEdz, E, zNow, zMax, relError);

  return std::min(value, VERYLARGEENERGY);
}

double Beniamino::dilationFactor(double E, double zNow, double zMax, double relError) const {
  assert(zMax >= zNow);
  if (generationEnergy(E, zNow, zMax, 1e-2) > 0.2 * VERYLARGEENERGY) return VERYLARGEJACOBIAN;
  auto dydz = [this, E, zNow](double z, double y) {
    auto dtdz = m_cosmology->dtdz(z);
    auto E_g = generationEnergy(E, zNow, z, 1e-6);
    auto E_prime = std::min(E_g * (1. + z) / (1. + zNow), VERYLARGEENERGY);
    auto dbdE = std::max(this->dbdE(E_prime), 0.);
    return y * (1. / (1. + z) + dtdz * pow3(1. + z) * dbdE);
  };
  auto value = utils::odeiv<double>(dydz, 1., zNow, zMax, relError);

  return std::min(value, VERYLARGEJACOBIAN);
}

// double Beniamino::computeFluxUnm(double E, double zMax, double relError) const {
//   const auto K = m_slope - 2;
//   const auto L_0 = m_sourceEmissivity;
//   const auto factor = SI::cLight / 4. / M_PI * K * L_0 * std::pow(E / E_0, -m_slope);
//   auto integrand = [&](double z) { return m_cosmology->dtdz(z) * std::pow(1. + z, -m_slope + 1.);
//   }; auto I = utils::QAGIntegration<double>(integrand, 0., zMax, INTSTEPS, relError);

//   return factor * I;
// }

double Beniamino::computeFlux(double E, double zObs, double zMax, double relError) const {
  const auto K = (m_slope - 2.) / pow2(m_minEnergy);
  const auto L_0 = m_sourceEmissivity;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0;
  auto integrand = [this, E, zObs](double z) {
    const auto E_g = generationEnergy(E, zObs, z, 1e-4);
    if (E_g > m_maxEnergy) return 0.;
    auto dEgdE = dilationFactor(E, zObs, z, 1e-4);
    auto inj = std::pow(E_g / m_minEnergy, -m_slope);
    if (m_sourceCutoff > 0.) inj *= std::exp(-E_g / m_sourceCutoff);
    auto sourceEvolution = std::pow(1. + z, m_sourceEvolution);
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz * sourceEvolution * inj * dEgdE;
  };
  auto I = utils::QAGIntegration<double>(integrand, zObs, zMax, INTSTEPS, relError);

  return factor * I;
}

}  // namespace solutions
}  // namespace simprop