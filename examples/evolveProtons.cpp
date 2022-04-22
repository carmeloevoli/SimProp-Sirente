#include <algorithm>

#include "simprop.h"

using namespace simprop;

auto IsActive = [](const Particle& p) {
  const double minPropagatingGamma = 1e8;
  return (p.getRedshift() > 1e-20 && p.getGamma() > minPropagatingGamma);
};

// auto IsPrimary = [](const Particle& p) { return p.IsPrimary(); };

template <typename T>
bool essentiallyEqual(T a, T b, T epsilon) {
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

class Evolutor {
 protected:
  const double deltaGammaCritical = 0.1;
  RandomNumberGenerator& m_rng;
  ParticleStack m_stack;
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::CMB> m_cmb;
  std::shared_ptr<photonfields::PhotonField> m_ebl;
  std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;
  std::shared_ptr<interactions::PhotoPionProduction> m_pppcmb;

 public:
  Evolutor(RandomNumberGenerator& rng) : m_rng(rng) {
    m_cosmology = std::make_shared<cosmo::Planck2018>();
  }

  void buildParticleStack(double z, double Gamma) {
    auto builder = SingleParticleBuilder(proton, 100);
    builder.setRedshift(z);
    builder.setGamma(Gamma);
    m_stack = builder.build(m_rng);
  }

  void buildPhotonFields() {
    m_cmb = std::make_shared<photonfields::CMB>();
    m_ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  }

  void buildContinuousLosses() {
    std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{m_cmb, m_ebl};
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::PairProductionLosses>(phFields),
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology)};
  }

  void buildStochasticInteractions() {
    auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
    m_pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, m_cmb);
  }

  double computeStochasticRedshiftInterval(PID pid, double Gamma, double zNow) {
    const auto dtdz = m_cosmology->dtdz(zNow);
    const auto lambda_s = std::fabs(1. / m_pppcmb->rate(pid, Gamma, zNow) / dtdz);
    // TODO why to put the fabs?
    return -lambda_s * std::log(1. - m_rng());
  }

  double computeDeltaGamma(PID pid, double Gamma, double zNow, double dz) {
    const auto dtdz = m_cosmology->dtdz(zNow);
    double dlnGammaNow = 0, dlnGammaHalf = 0, dlnGammaNext = 0;
    for (auto losses : m_continuousLosses) {
      dlnGammaNow += losses->dlnGamma_dt(pid, Gamma, zNow);
      auto halfRedshift = zNow - 0.5 * dz;
      dlnGammaHalf += losses->dlnGamma_dt(pid, Gamma, halfRedshift);
      auto nextRedhisft = zNow - dz;
      dlnGammaNext += losses->dlnGamma_dt(pid, Gamma, nextRedhisft);
    }
    return dz / 6. * dtdz * (dlnGammaNow + 4. * dlnGammaHalf + dlnGammaNext);
  }

  double computeLossesRedshiftInterval(PID pid, double Gamma, double zNow) {
    double dz = zNow;
    double deltaGamma = computeDeltaGamma(pid, Gamma, zNow, zNow);
    if (deltaGamma > deltaGammaCritical) {
      dz = utils::rootFinder<double>(
          [&](double x) { return computeDeltaGamma(pid, Gamma, zNow, x) - deltaGammaCritical; }, 0.,
          zNow, 100, 1e-5);
    }
    return dz;
  }

  void run(std::string filename) {
    ParticleStack::iterator it = m_stack.begin();
    utils::OutputFile out(filename.c_str());
    while (it != m_stack.end()) {
      const auto pid = it->getPid();
      const auto nowRedshift = it->getRedshift();
      const auto Gamma = it->getGamma();

      const auto dz_s = computeStochasticRedshiftInterval(pid, Gamma, nowRedshift);
      assert(dz_s > 0.);

      const auto dz_c = computeLossesRedshiftInterval(pid, Gamma, nowRedshift);
      assert(dz_c > 0. && dz_c <= nowRedshift);

      if (dz_s > dz_c || dz_s > nowRedshift) {
        const auto dz = dz_c;
        it->getNow().z -= dz;
        auto deltaGamma = computeDeltaGamma(pid, Gamma, nowRedshift, dz);
        it->getNow().Gamma *= (1. - deltaGamma);
#ifdef PRINTALL
        out << *it << " " << 0 << "\n";
#endif
      } else {
        const auto dz = dz_s;
        it->getNow().z -= dz;
        auto state = m_pppcmb->finalState(*it, nowRedshift - dz, m_rng);
        it->getNow().Gamma = state.at(0).getGamma();
#ifdef PRINTALL
        out << *it << " " << 1 << "\n";
#endif
      }

#ifdef PRINTALL
      std::cout << *it << " " << dz_s << " " << dz_c << "\n";
#endif

      it = std::find_if(m_stack.begin(), m_stack.end(), IsActive);
    }
  }
  virtual ~Evolutor() = default;
};

int main() {
  try {
    utils::startup_information();
    RandomNumberGenerator rng = utils::RNG<double>(66);
    {
      utils::Timer timer("timer for Gamma = 1e12");
      Evolutor evolutor(rng);
      evolutor.buildParticleStack(1., 1e12);
      evolutor.buildPhotonFields();
      evolutor.buildContinuousLosses();
      evolutor.buildStochasticInteractions();
      evolutor.run("test_proton_evolution_1e12.txt");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
