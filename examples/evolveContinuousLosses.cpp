#include <algorithm>

#include "simprop.h"

using namespace simprop;

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

  void buildParticleStack(Redshift z, LorentzFactor Gamma, size_t N = 1) {
    auto builder = SingleParticleBuilder(proton, {Gamma.get(), z.get()}, N);
    m_stack = builder.build();
  }

  void buildCosmologicalParticleStack(double slope, double evolutionIndex, Redshift zMax,
                                      size_t N = 1) {
    double maxEnergy = std::pow(10., 23) * SI::eV;
    double minEnergy = std::pow(10., 17) * SI::eV;
    Range GammaRange = {minEnergy / SI::protonMassC2, maxEnergy / SI::protonMassC2};
    Range zRange = {0., zMax.get()};
    auto builder =
        SourceEvolutionBuilder(proton, {GammaRange, zRange, slope, evolutionIndex}, m_cosmology, N);
    m_stack = builder.build(m_rng);
  }

  void buildPhotonFields() { m_cmb = std::make_shared<photonfields::CMB>(); }

  void buildContinuousLosses(bool doPionLosses = false) {
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::BGG2002ContinuousLosses>(),
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology)};
    if (doPionLosses)
      m_continuousLosses.push_back(std::make_shared<losses::PhotoPionContinuousLosses>());
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
          100, 1e-5);
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
      // const double dz_c = 1e-3 * (1. + nowRedshift);
      const auto dz_c = computeLossesRedshiftInterval(*it);
      assert(dz_c > 0. && dz_c <= nowRedshift);
      const auto deltaGamma = computeDeltaGamma(*it, dz_c);
      it->getNow() = {nowRedshift - dz_c, Gamma * (1. - deltaGamma)};
      it = std::find_if(it, m_stack.end(), IsActive);
      counter++;
    }
  }  // run()

  double getObservedEnergy() {  // TODO remove this
    assert(m_stack.size() == 1 && m_stack[0].getRedshift() < 1e-20);
    return m_stack[0].getGamma() * SI::protonMassC2;
  }

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

  virtual ~Evolutor() = default;
};

void testSingleProtonEvolution() {
  RandomNumberGenerator rng = utils::RNG<double>(Seed(96));
  const double zMax = 3.0;
  utils::OutputFile out("test_traj_z3.0.txt");
  for (double E = 1e17 * SI::eV; E < 1e23 * SI::eV; E *= 1.25) {
    Evolutor evolutor(rng);
    evolutor.buildParticleStack(Redshift(zMax), LorentzFactor(E / SI::protonMassC2), 1);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses(false);
    evolutor.run("test_proton_evolution_1_1e17.txt");
    out << std::scientific << E / SI::eV << " " << evolutor.getObservedEnergy() / SI::eV << "\n";
  }
}

void evolvePopulation() {
  RandomNumberGenerator rng = utils::RNG<double>(Seed(69));
  {
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.7, 0., Redshift(0.05), 1e6);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_spectrum_z0.05_m0_smart_N1e6.txt");
  }
  {
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.7, 0., Redshift(0.5), 1e6);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_spectrum_z0.5_m0_smart_N1e6.txt");
  }
  {
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.7, 0., Redshift(1.0), 1e6);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_spectrum_z1.0_m0_smart_N1e6.txt");
  }
  {
    utils::Timer timer("timer for z=3");
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.7, 0., Redshift(3.0), 1e6);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_spectrum_z3.0_m0_smart_N1e6.txt");
  }
}

int main() {
  try {
    utils::startup_information();
    // testSingleProtonEvolution();
    evolvePopulation();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
