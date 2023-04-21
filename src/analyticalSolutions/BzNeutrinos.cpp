#include "simprop/analyticalSolutions/BzNeutrinos.h"

#include <limits>

#include "simprop/analyticalSolutions/beniamino.h"
#include "simprop/utils/logging.h"

namespace simprop {
namespace solutions {

BzNeutrinos::BzNeutrinos() {
  m_cosmology = std::make_shared<cosmo::Cosmology>();
  m_nuSpec = std::make_shared<KelnerAharonian2008::NeutrinoProductionSpectrum>();
  m_cmb = std::make_shared<photonfields::CMB>();

  using std::exp;
  using std::log;
  {
    solutions::Beniamino b(true);
    const auto zMax = 6.;
    auto f = [&b, zMax](double logEp, double z) -> double {
      return log(b.computeFlux(exp(logEp), z, zMax, 1e-2));
    };
    m_Jp.cacheTable(f, {log(1e16 * SI::eV), log(1e23 * SI::eV)}, {0., zMax});
  }
}

double BzNeutrinos::I_deps(double Enu, double Ep, double z) const {
  const auto epsMin = m_cmb->getMinPhotonEnergy();
  const auto epsMax = m_cmb->getMaxPhotonEnergy();
  auto integrand = [this, Ep, Enu, z](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = m_cmb->density(epsilon, z) * m_nuSpec->Phi(eta, Enu / Ep);
    return epsilon * value;
  };
  auto I = utils::simpsonIntegration<double>(integrand, std::log(epsMin), std::log(epsMax), 300);
  return I;
}

double BzNeutrinos::I_dEp(double Enu, double z) const {
  auto integrand = [this, Enu, z](double lnEp) {
    const auto Ep = std::exp(lnEp);
    auto value = m_Jp.get(lnEp, z);
    value *= I_deps(Enu, Ep, z);
    return value;
  };

  return utils::QAGIntegration<double>(integrand, std::log(Enu), std::log(1e5 * Enu), 1000, 1e-2);
}

double BzNeutrinos::computeNeutrinoFlux(double EnuObs, double zMax, double relError) const {
  const auto factor = 1.;
  auto integrand = [this, EnuObs](double z) {
    auto dtdz = m_cosmology->dtdz(z);
    return dtdz / pow2(1. + z) * I_dEp(EnuObs * (1. + z), z);
  };
  auto I = utils::simpsonIntegration<double>(integrand, 0., zMax, 100);

  return factor * I;
}

double BzNeutrinos::getProtonFlux(double Ep, double z) const {
  auto lnEp = std::log(Ep);
  return std::exp(m_Jp.get(lnEp, z));
}

}  // namespace solutions
}  // namespace simprop