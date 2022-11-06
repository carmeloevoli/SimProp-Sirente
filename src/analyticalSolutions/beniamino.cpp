#include "simprop/analyticalSolutions/beniamino.h"

#include <limits>

namespace simprop {
namespace solutions {

#define INTSTEPS 2000

Beniamino::Beniamino() {
  m_cosmology = std::make_shared<cosmo::Planck2018>();
  m_pp = std::make_shared<losses::BGG2002ContinuousLosses>();
  m_pion = std::make_shared<losses::PhotoPionContinuousLosses>();
}

double Beniamino::generationEnergy(double E, double zMax, double relError) const {
  auto dEdz = [&](double z, double E_g) {
    auto Gamma = E_g / SI::protonMassC2 * (1. + z);
    auto beta = m_pp->beta(proton, Gamma) + ((m_doPhotoPion) ? m_pion->beta(proton, Gamma) : 0.);
    auto dtdz = m_cosmology->dtdz(z);
    return E_g * (1. / (1. + z) + dtdz * pow3(1. + z) * beta);
  };
  auto value = utils::odeiv<double>(dEdz, E, 0., zMax, relError);
  return std::min(value, std::numeric_limits<double>::max());
}

double Beniamino::dbdE(double E) const {
  auto dbeta = utils::deriv5pt<double>(
      [&](double x) {
        auto Gamma = x / SI::protonMassC2;
        auto beta = m_pp->beta(proton, Gamma);
        beta += ((m_doPhotoPion) ? m_pion->beta(proton, Gamma) : 0.);
        return x * beta;
      },
      E, 1e-2 * E);
  return dbeta;
}

double Beniamino::dilationFactor(double E, double zMax, double relError) const {
  auto integrand = [&](double z) {
    auto E_g = generationEnergy(E, z, 1e-6);
    auto E_prime = E_g * (1. + z);
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz * pow3(1. + z) * dbdE(E_prime);
  };
  auto I = utils::simpsonIntegration<double>(integrand, 0., zMax, INTSTEPS);
  //  auto I = utils::QAGSIntegration<double>(integrand, 0., zMax, INTSTEPS, 1e-3);
  auto value = (1. + zMax) * std::exp(I);
  return (value < 1.) ? 1e10 : value;
}

double Beniamino::computeFluxUnm(double E) const {
  const auto K = m_slope - 2;
  const auto L_0 = m_sourceEmissivity;
  const auto E_0 = 1e17 * SI::eV;
  const auto zMax = m_sourceMaxRedshift;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0 * std::pow(E / E_0, -m_slope);
  auto integrand = [&](double z) { return m_cosmology->dtdz(z) * std::pow(1. + z, -m_slope + 1.); };
  // auto I = utils::simpsonIntegration<double>(integrand, 0., zMax, INTSTEPS);
  auto I = utils::QAGIntegration<double>(integrand, 0., zMax, INTSTEPS, 1e-4);
  return factor * I;
}

double Beniamino::computeFlux(double E) const {
  const auto K = m_slope - 2;
  const auto L_0 = m_sourceEmissivity;
  const auto m = m_sourceEvolution;
  const auto E_0 = 1e17 * SI::eV;
  const auto zMax = m_sourceMaxRedshift;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0;
  auto integrand = [&](double z) {
    auto E_g = generationEnergy(E, z, 1e-5);
    if (E_g > m_maxEnergy) return 0.;
    auto dEgdE = dilationFactor(E, z, 1e-5);
    auto inj = std::pow(E_g / E_0, -m_slope);
    if (m_sourceCutoff > 0.) inj *= std::exp(-E_g / m_sourceCutoff);
    auto sourceEvolution = std::pow(1. + z, m);
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz * sourceEvolution * inj * dEgdE;
  };
  // auto I = utils::simpsonIntegration<double>(integrand, 0., zMax, INTSTEPS);
  auto I = utils::QAGIntegration<double>(integrand, 0., zMax, INTSTEPS, 1e-3);
  return factor * I;
}

}  // namespace solutions
}  // namespace simprop