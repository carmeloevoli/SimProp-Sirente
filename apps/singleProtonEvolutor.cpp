#include "simprop.h"

using namespace simprop;

#define VERYLARGEENERGY (1e25 * SI::eV)

auto IsActive = [](const Particle& p) {
  constexpr double minPropagatingGamma = 1e6;
  return (p.isNucleus() && p.isActive() && p.getRedshift() > 1e-20 &&
          p.getGamma() > minPropagatingGamma);
};

class Evolutor {
 protected:
  const double deltaGammaCritical = 0.03;
  RandomNumberGenerator& m_rng;
  ParticleStack m_stack;
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_cmb;
  std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;

 public:
  explicit Evolutor(RandomNumberGenerator& rng) : m_rng(rng) {}

  inline const ParticleStack& getStack() const { return m_stack; }

  void buildParticleStack(double z, double Gamma, size_t N = 1) {
    auto builder = SingleParticleBuilder(proton, {Gamma, z}, N);
    m_stack = builder.build();
  }

  void clear() { m_stack.clear(); }

  void buildCosmology() { m_cosmology = std::make_shared<cosmo::Planck2018>(); }

  void buildPhotonFields() { m_cmb = std::make_shared<photonfields::CMB>(); }

  void buildContinuousLosses(bool doPionLosses = false) {
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology),
        std::make_shared<losses::PairProductionLosses>(m_cmb),
    };
    if (doPionLosses)
      m_continuousLosses.push_back(std::make_shared<losses::PhotoPionContinuousLosses>(m_cmb));
  }

  double computeDeltaGamma(const Particle& particle, double deltaRedshift) const {
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

  double computeLossesRedshiftInterval(const Particle& particle) const {
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

  void run() {
    size_t counter = 0;
    auto it = m_stack.begin();
    while (it != m_stack.end()) {
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

  double getObservedEnergy() {
    assert(m_stack.size() == 1);  // && m_stack[0].getRedshift() < 1e-20);
    return m_stack[0].getGamma() * SI::protonMassC2;
  }
};

void testSingleParticleEvolution(double zMax, std::string filename) {
  RandomNumberGenerator rng = utils::RNG<double>(69);
  utils::OutputFile out(filename);
  Evolutor evolutor(rng);
  evolutor.buildCosmology();
  evolutor.buildPhotonFields();
  evolutor.buildContinuousLosses(true);
  const auto eFactor = std::pow(10., 1. / 4.);
  for (double E = 1e16 * SI::eV; E < 1e23 * SI::eV; E *= eFactor) {
    evolutor.buildParticleStack(zMax, E / SI::protonMassC2, 1);
    evolutor.run();
    out << std::scientific << E / SI::eV << " " << evolutor.getObservedEnergy() / SI::eV << "\n";
    evolutor.clear();
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    testSingleParticleEvolution(5.0, "test_traj_z5.0.txt");
    testSingleParticleEvolution(3.0, "test_traj_z3.0.txt");
    testSingleParticleEvolution(2.0, "test_traj_z2.0.txt");
    testSingleParticleEvolution(1.0, "test_traj_z1.0.txt");
    testSingleParticleEvolution(0.5, "test_traj_z0.5.txt");
    testSingleParticleEvolution(0.05, "test_traj_z0.05.txt");
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}