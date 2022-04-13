#include <algorithm>

#include "simprop.h"

using namespace simprop;

bool isActive(Particle p) {
  const double minPropagatingGamma = 1e8;
  return (p.getRedshift() > 1e-20 && p.getGamma() > minPropagatingGamma);
}

template <typename T>
bool essentiallyEqual(T a, T b, T epsilon) {
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
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
    auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);

    ParticleStack::iterator it = stack.begin();
    utils::OutputFile out("test_proton_evolution_1e10.txt");

    while (it != stack.end()) {
      const auto nowRedshift = it->getNow().z;
      const auto Gamma = it->getNow().Gamma;
      const auto dtdz = cosmology->dtdz(nowRedshift);

      auto lambda_s = std::fabs(1. / pppcmb->rate(it->getPid(), Gamma, nowRedshift) /
                                dtdz);  // TODO check this fabs!
      auto dz_s = -lambda_s * std::log(1. - rng());

      // std::cout << *it << " dz_s : " << dz_s << " lambda_s : " << lambda_s << "\n";
      assert(dz_s > 0.);

      auto dz_c = nowRedshift;
      double deltaGamma = 100.0;
      {
        double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
        for (auto losses : continuousLosses) {
          dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
          auto halfRedshift = nowRedshift - 0.5 * dz_c;
          dlnGammaHalf += losses->dlnGamma_dt(it->getPid(), Gamma, halfRedshift);
          auto nextRedhisft = nowRedshift - dz_c;
          dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedhisft);
        }
        deltaGamma = dz_c / 6. * dtdz * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext);
      }
      while (deltaGamma > 0.1) {  // TODO better a bisection method?
        dz_c *= 0.9;
        double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
        for (auto losses : continuousLosses) {
          dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
          auto halfRedshift = nowRedshift - 0.5 * dz_c;
          dlnGammaHalf += losses->dlnGamma_dt(it->getPid(), Gamma, halfRedshift);
          auto nextRedhisft = nowRedshift - dz_c;
          dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedhisft);
        }
        deltaGamma = dz_c / 6. * dtdz * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext);
      }

      // std::cout << *it << " dz_c : " << dz_c << " z_now : " << nowRedshift << "\n";
      assert(dz_c > 0. && dz_c <= nowRedshift);

      if (dz_s > dz_c || dz_s > nowRedshift) {
        const auto dz = dz_c;
        const auto dt = dtdz * dz;
        double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
        for (auto losses : continuousLosses) {
          dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
          auto halfRedshift = nowRedshift - 0.5 * dz_c;
          dlnGammaHalf += losses->dlnGamma_dt(it->getPid(), Gamma, halfRedshift);
          auto nextRedhisft = nowRedshift - dz_c;
          dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedhisft);
        }
        it->getNow().Gamma *= 1. - dt / 6. * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext);
        it->getNow().z -= dz;
        out << *it << " " << 0 << "\n";
      } else {
        const auto dz = dz_s;
        const auto dt = dtdz * dz;
        double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
        for (auto losses : continuousLosses) {
          dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
          auto halfRedshift = nowRedshift - 0.5 * dz_c;
          dlnGammaHalf += losses->dlnGamma_dt(it->getPid(), Gamma, halfRedshift);
          auto nextRedhisft = nowRedshift - dz_c;
          dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedhisft);
        }
        it->getNow().Gamma *= 1. - dt / 6. * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext);
        it->getNow().z -= dz;
        auto state = pppcmb->finalState(*it, nowRedshift - dz, rng);
        it->getNow().Gamma = state.at(0).getGamma();
        out << *it << " " << 1 << "\n";
      }

      //   const double dz = 0.01;
      //   const auto nextRedshift = std::max(nowRedshift - dz, 0.);
      //   it->getNow().z = nextRedshift;

      //   //   if (DeltazInteraction < dz) {
      //   //     auto state = pppcmb->finalState(*it, nowRedshift - DeltazInteraction, rng);
      //   //     it->getNow().Gamma = state.at(0).getGamma();
      //   //   }

      //   Gamma = it->getNow().Gamma;
      //   double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
      //   for (auto losses : continuousLosses) {
      //     dlnGammaNow += losses->dlnGamma_dt(it->getPid(), Gamma, nowRedshift);
      //     dlnGammaHalf +=
      //         losses->dlnGamma_dt(it->getPid(), Gamma, 0.5 * (nowRedshift + nextRedshift));
      //     dlnGammaNext += losses->dlnGamma_dt(it->getPid(), Gamma, nextRedshift);
      //   }

      //   const auto dt = cosmology->dtdz(nowRedshift) * dz;
      //   it->getNow().Gamma =
      //       Gamma * (1. - dt / 6. * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext));

      //   //   if (dtInteractionPoint < dt)
      //   //     std::cout << "interaction take place"
      //   //               << "\n";

      //   //    = std::accumulate(continuousLosses.begin()), continuousLosses.end(),
      //   //                                0, [Gamma](std::shared_ptr<losses::ContinuousLosses>
      //   l0,
      //   //                                std::shared_ptr<losses::ContinuousLosses> l1) {
      //   //     return l0 + l1;});

      //   //   const auto nowEnergy = it->getNow().E;
      //   //   it->getNow().E =
      //   //       m_continuousLosses->evolve(nowEnergy, nowRedshift, nextRedshift, it->getPid());
      //   out << *it << "\n";
      it = std::find_if(stack.begin(), stack.end(), isActive);
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
