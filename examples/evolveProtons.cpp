#include <algorithm>

#include "simprop.h"

using namespace simprop;

bool isActive(Particle p) {
  const double minPropagatingGamma = 1e8;
  return (p.getRedshift() > 1e-20 && p.getGamma() > minPropagatingGamma);
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    RandomNumberGenerator rng = utils::RNG<double>(66);
    auto stack = ParticleStack(proton, 1, rng);
    stack.buildSingleParticleStack(1., 1e10);
    auto cosmology = std::make_shared<cosmo::Planck2018>();

    std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{
        std::make_shared<photonfields::CMB>(),
        std::make_shared<photonfields::Dominguez2011PhotonField>()};

    std::vector<std::shared_ptr<losses::ContinuousLosses> > continuousLosses{
        std::make_shared<losses::PairProductionLosses>(phFields),
        std::make_shared<losses::AdiabaticContinuousLosses>(cosmology)};

    const auto cmb = std::make_shared<photonfields::CMB>();
    auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
    auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(
        sigma, cmb);  // TODO pass phFields here!

    ParticleStack::iterator it = stack.begin();
    while (it != stack.end()) {
      const double dz = 0.001;
      const auto nowRedshift = it->getNow().z;
      const auto nextRedshift = std::max(nowRedshift - dz, 0.);
      it->getNow().z = nextRedshift;

      const auto Gamma = it->getNow().Gamma;
      double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
      for (auto losses : continuousLosses) {
        dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
        dlnGammaHalf +=
            losses->dlnGamma_dt(it->getPid(), Gamma, 0.5 * (nowRedshift + nextRedshift));
        dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedshift);
      }

      const auto dt = cosmology->dtdz(nowRedshift);
      it->getNow().Gamma =
          Gamma * (1. - dt / 6. * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext));

      //    = std::accumulate(continuousLosses.begin()), continuousLosses.end(),
      //                                0, [Gamma](std::shared_ptr<losses::ContinuousLosses> l0,
      //                                std::shared_ptr<losses::ContinuousLosses> l1) {
      //     return l0 + l1;});

      //   const auto nowEnergy = it->getNow().E;
      //   it->getNow().E =
      //       m_continuousLosses->evolve(nowEnergy, nowRedshift, nextRedshift, it->getPid());
      std::cout << *it << "\n";
      it = std::find_if(stack.begin(), stack.end(), isActive);
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
