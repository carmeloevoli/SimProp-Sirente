#include "simprop.h"

using namespace simprop;

#define VERYLARGEENERGY (1e25 * SI::eV)

auto IsActive = [](const Particle& p) {
  constexpr double minPropagatingGamma = 1e7;
  return (p.isNucleus() && p.isActive() && p.getRedshift() > 1e-20 &&
          p.getGamma() > minPropagatingGamma);
};

class Evolutor {
 protected:
  const double deltaGammaCritical = 0.1;
  RandomNumberGenerator& m_rng;
  ParticleStack m_stack;
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_cmb;
  std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;

 public:
  Evolutor(RandomNumberGenerator& rng) : m_rng(rng) {
    m_cosmology = std::make_shared<cosmo::Planck2018>();
  }

  inline const ParticleStack& getStack() const { return m_stack; }

  void buildParticleStack(double z, double Gamma, size_t N = 1) {
    auto builder = SingleParticleBuilder(proton, {Gamma, z}, N);
    m_stack = builder.build();
  }

  void buildCosmologicalParticleStack(double slope, double evolutionIndex, double zMax,
                                      size_t N = 1) {
    double maxEnergy = std::pow(10., 22) * SI::eV;
    double minEnergy = std::pow(10., 17) * SI::eV;
    Range GammaRange = {minEnergy / SI::protonMassC2, maxEnergy / SI::protonMassC2};
    Range zRange = {0., zMax};
    auto builder =
        SourceEvolutionBuilder(proton, {GammaRange, zRange, slope, evolutionIndex}, m_cosmology, N);
    m_stack = builder.build(m_rng);
  }

  void buildPhotonFields() { m_cmb = std::make_shared<photonfields::CMB>(); }

  void buildContinuousLosses(bool doPionLosses = false) {
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology),
        std::make_shared<losses::PairProductionLosses>(
            losses::PairProductionLosses(m_cmb).doCaching()),
    };
    if (doPionLosses)
      m_continuousLosses.push_back(std::make_shared<losses::PhotoPionContinuousLosses>(m_cmb));
  }

  double computeDeltaGamma(const Particle& particle, double deltaRedshift) {
    const auto pid = particle.getPid();
    const auto Gamma = particle.getGamma();
    const auto zNow = particle.getRedshift();
    const auto zHalf = zNow - 0.5 * deltaRedshift;
    const auto zNext = zNow - deltaRedshift;
    double betaNow = 0, betaHalf = 0, betaNext = 0;
    for (auto losses : m_continuousLosses) {
      betaNow += losses->beta(pid, Gamma, zNow);
      betaHalf += losses->beta(pid, Gamma, zHalf);
      betaNext += losses->beta(pid, Gamma, zNext);
    }
    double value = deltaRedshift / 6.;
    value *= betaNow * m_cosmology->dtdz(zNow) + 4. * betaHalf * m_cosmology->dtdz(zHalf) +
             betaNext * m_cosmology->dtdz(zNext);
    return 1.0 - std::exp(-value);
  }

  double computeLossesRedshiftInterval(const Particle& particle) {
    const auto zNow = particle.getRedshift();
    double dz = zNow;

    double deltaGamma = computeDeltaGamma(particle, zNow);
    if (deltaGamma > deltaGammaCritical) {
      dz = utils::rootFinder<double>(
          [&](double x) { return computeDeltaGamma(particle, x) - deltaGammaCritical; }, 0., zNow,
          100, 1e-3);
    }
    return dz;
  }

  void run(std::string filename) {
    utils::OutputFile out(filename.c_str());
    // size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
    size_t iniSize = m_stack.size();
    const double initEnergy = sumEnergy();
    size_t counter = 0;
    auto it = m_stack.begin();
    while (it != m_stack.end()) {
      if (counter % iniSize == 0) {
        LOGD << counter / iniSize << "\t" << sumEnergy() / initEnergy << "\t" << m_stack.size();
      }
      const auto nowRedshift = it->getRedshift();
      const auto Gamma = it->getGamma();
      const auto dz_c = computeLossesRedshiftInterval(*it);
      assert(dz_c > 0. && dz_c <= nowRedshift);
      const auto deltaGamma = computeDeltaGamma(*it, dz_c);
      it->getNow() = {nowRedshift - dz_c, Gamma * (1. - deltaGamma)};
      it = std::find_if(it, m_stack.end(), IsActive);
      counter++;
    }
  }  // run()

  // double getObservedEnergy() {  // TODO remove this
  //   assert(m_stack.size() == 1 && m_stack[0].getRedshift() < 1e-20);
  //   return m_stack[0].getGamma() * SI::protonMassC2;
  // }

  double sumEnergy() {
    double value = 0;
    std::for_each(m_stack.begin(), m_stack.end(), [&](const auto& particle) {
      value += (IsActive(particle)) ? particle.getEnergy() : 0.;
    });
    return value;
  }

  void dump(const std::string& filename) {
    utils::OutputFile out(filename.c_str());
    for (const auto& particle : m_stack) {
      if (particle.getPid() == proton && particle.getRedshift() < 1e-20) out << particle << "\n";
    }
  }
};

// std::shared_ptr<interactions::PhotoPionProduction> m_pppcmb;

// void buildStochasticInteractions() {
//   auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
//   m_pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, m_cmb);
// }

// double computeStochasticRedshiftInterval(const Particle& particle) {
//   const auto pid = particle.getPid();
//   const auto Gamma = particle.getGamma();
//   const auto zNow = particle.getRedshift();
//   const auto dtdz = m_cosmology->dtdz(zNow);
//   const auto dt = std::fabs(1. / m_pppcmb->rate(pid, Gamma, zNow));
//   // TODO why to put the fabs?
//   return -dt / dtdz * std::log(1. - m_rng());
// }

// void run(std::string filename) {
//   assert(m_pppcmb != nullptr);
//   auto it = m_stack.begin();
//   size_t iniSize = m_stack.size();
//   size_t counter = 0;
//   const double initEnergy = sumEnergy();
//   while (it != m_stack.end()) {
//     if (counter % iniSize == 0) {
//       LOGD << counter / iniSize << "\t" << sumEnergy() / initEnergy << "\t" << m_stack.size();
//     }

//     const auto distance = it - m_stack.begin();
//     const auto nowRedshift = it->getRedshift();
//     const auto dz_s = computeStochasticRedshiftInterval(*it);
//     const auto dz_c = computeLossesRedshiftInterval(*it);
//     assert(dz_s > 0. && dz_c > 0. && dz_c <= nowRedshift);

//     if (dz_s > dz_c || dz_s > nowRedshift) {
//       const auto Gamma = it->getGamma();
//       const auto dz = dz_c;
//       const auto deltaGamma = computeDeltaGamma(*it, dz);
//       it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
//     } else {
//       it->deactivate();
//       const auto dz = dz_s;
//       auto finalState = m_pppcmb->finalState(*it, nowRedshift - dz, m_rng);
//       m_stack.insert(m_stack.end(), finalState.begin(), finalState.end());
//     }
//     it = std::find_if(m_stack.begin() + distance, m_stack.end(), IsActive);
//     counter++;
//   }
// }  // run()

void evolvePopulation() {
  RandomNumberGenerator rng = utils::RNG<double>(69);
  {
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.6, 0., 1.0, 1e5);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses(true);
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_spectrum_z1.0_m0_N1e6.txt");
  }
}

void testSpectrum() {
  solutions::Beniamino b(true);
  auto E = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 16 * 5);
  {
    utils::OutputFile out("test_analytical_spectrum.txt");
    const double zMax = 1.;
    const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
    for (const auto& E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i, 0., zMax, 1e-2) / units << "\t";
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    testSpectrum();
    // evolvePopulation();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}