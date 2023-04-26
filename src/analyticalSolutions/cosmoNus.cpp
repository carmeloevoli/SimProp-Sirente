#include "simprop/analyticalSolutions/cosmoNus.h"

#include <limits>

#include "simprop/utils/logging.h"

namespace simprop {
namespace solutions {

CosmoNeutrinos::CosmoNeutrinos(double injSlope, double sourceEvolution, double sourceCutoff) {
  m_cosmology = std::make_shared<cosmo::Cosmology>();
  m_nuSpec = std::make_shared<KelnerAharonian2008::NeutrinoProductionSpectrum>();
  m_cmb = std::make_shared<photonfields::CMB>();
  {
    using std::exp;
    using std::log;
    auto b = solutions::Beniamino(injSlope, sourceEvolution, sourceCutoff, true).doCaching();
    const auto zMax = 5.;
    auto f = [&b, zMax](double logEp, double z) -> double {
      auto value = b.computeFlux(exp(logEp), z, zMax, 1e-2);
      return log(std::max(value, 1e-30));
    };
    m_Jp.cacheTable(f, {log(1e16 * SI::eV), log(1e23 * SI::eV)}, {0., zMax});
  }
}

double CosmoNeutrinos::I_deps(double Enu, double Ep, double z, size_t N) const {
  const auto epsMin = m_cmb->getMinPhotonEnergy();
  const auto epsMax = m_cmb->getMaxPhotonEnergy();
  auto integrand = [this, Ep, Enu, z](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = m_cmb->density(epsilon, z) * m_nuSpec->Phi(eta, Enu / Ep);
    return epsilon * value;
  };
  auto a = std::log(epsMin);
  auto b = std::log(epsMax);
  auto I = utils::RombergIntegration<double>(integrand, a, b, N, 1e-2);

  return I;
}

double CosmoNeutrinos::I_dEp(double Enu, double z, size_t N) const {
  auto integrand = [this, Enu, z](double lnEp) {
    const auto Ep = std::exp(lnEp);
    auto value = std::exp(m_Jp.get(lnEp, z));
    value *= I_deps(Enu, Ep, z);
    return value;
  };
  auto a = std::log(Enu);
  auto b = std::log(std::min(1e4 * Enu, 1e23 * SI::eV));
  auto I = utils::RombergIntegration<double>(integrand, a, b, N, 1e-2);

  return I;
}

double CosmoNeutrinos::computeNeutrinoFlux(double EnuObs, double zMax, size_t N) const {
  const auto factor = 1.;
  auto integrand = [this, EnuObs](double z) {
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz / pow2(1. + z) * I_dEp(EnuObs * (1. + z), z);
  };
  auto I = utils::RombergIntegration<double>(integrand, 0., zMax, N, 1e-2);

  return factor * I;
}

double CosmoNeutrinos::getProtonFlux(double Ep, double z) const {
  auto lnEp = std::log(Ep);
  return std::exp(m_Jp.get(lnEp, z));
}

}  // namespace solutions
}  // namespace simprop