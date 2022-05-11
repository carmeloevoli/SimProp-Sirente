#include <algorithm>

#include "simprop.h"

using namespace simprop;

auto IsActive = [](const Particle& p) {
  constexpr double minPropagatingGamma = 1e7;
  return (p.IsNucleus() && p.getRedshift() > 1e-20 && p.getGamma() > minPropagatingGamma);
};

class Evolutor {
 protected:
  const double deltaGammaCritical = 0.1;
  RandomNumberGenerator& m_rng;
  ParticleStack m_stack;
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_cmb;
  // std::shared_ptr<photonfields::PhotonField> m_ebl;
  std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;
  std::shared_ptr<interactions::PhotoPionProduction> m_pppcmb;

 public:
  Evolutor(RandomNumberGenerator& rng) : m_rng(rng) {
    m_cosmology = std::make_shared<cosmo::Planck2018>();
  }

  void buildParticleStack(Redshift z, LorentzFactor Gamma, size_t N = 1) {
    auto builder = SingleParticleBuilder(proton, {Gamma.get(), z.get()}, N);
    m_stack = builder.build();
  }

  void buildCosmologicalParticleStack(double slope, double evolutionIndex, double zMax,
                                      size_t N = 1) {
    double maxEnergy = std::pow(10., 23) * SI::eV;
    double minEnergy = std::pow(10., 17) * SI::eV;
    Range GammaRange = {minEnergy / SI::protonMassC2, maxEnergy / SI::protonMassC2};
    Range zRange = {0., zMax};
    auto builder =
        SourceEvolutionBuilder(proton, {GammaRange, zRange, slope, evolutionIndex}, m_cosmology, N);
    m_stack = builder.build(m_rng);
  }

  void buildPhotonFields() {
    m_cmb = std::make_shared<photonfields::CMB>();
    // m_ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  }

  void buildContinuousLosses() {
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::BGG2002ContinuousLosses>(),
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology)};
  }

  void buildContinuousLossesIncludingPion() {
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology),
        std::make_shared<losses::BGG2002ContinuousLosses>(),
        std::make_shared<losses::PhotoPionContinuousLosses>()};
  }

  void buildStochasticInteractions() {
    auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
    m_pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, m_cmb);
  }

  double computeStochasticRedshiftInterval(const Particle& particle) {
    const auto pid = particle.getPid();
    const auto Gamma = particle.getGamma();
    const auto zNow = particle.getRedshift();
    const auto dtdz = m_cosmology->dtdz(zNow);
    const auto dt = std::fabs(1. / m_pppcmb->rate(pid, Gamma, zNow));
    // TODO why to put the fabs?
    return -dt / dtdz * std::log(1. - m_rng());
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

  double accumulateEnergy() {
    double value;
    std::for_each(m_stack.begin(), m_stack.end(), [&](const auto& particle) {
      value += (IsActive(particle)) ? particle.getEnergy() : 0.;
    });
    return value;
  }

  void run(std::string filename) {
    assert(m_pppcmb != nullptr);
    utils::OutputFile out(filename.c_str());
    size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
    // auto it = m_stack.begin();
    size_t counter = 0;
    // while (it != m_stack.end())
    while (nActive > 0) {
      if (counter % 1000 == 0) {
        std::cout << counter / 1000 << " " << accumulateEnergy() / SI::eV << "\n";
      }
      const auto it = std::find_if(m_stack.begin(), m_stack.end(), IsActive);
      // out << *it << " " << 0 << "\n";
      const auto nowRedshift = it->getRedshift();

      const auto dz_s = computeStochasticRedshiftInterval(*it);
      assert(dz_s > 0.);

      const auto dz_c = computeLossesRedshiftInterval(*it);
      assert(dz_c > 0. && dz_c <= nowRedshift);

      if (dz_s > dz_c || dz_s > nowRedshift) {
        const auto Gamma = it->getGamma();
        const auto dz = dz_c;
        const auto deltaGamma = computeDeltaGamma(*it, dz);
        it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
        // out << *it << " " << 0 << "\n";
      } else {
        // out << *it << " " << 1 << "\n";
        const auto dz = dz_s;
        auto finalState = m_pppcmb->finalState(*it, nowRedshift - dz, m_rng);
        m_stack.erase(it);
        m_stack.insert(m_stack.begin(), finalState.begin(), finalState.end());
      }

      nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
      counter++;
    }
  }  // run()

  void runContinuous(std::string filename) {
    utils::OutputFile out(filename.c_str());
    size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
    size_t counter = 0;
    auto it = m_stack.begin();
    while (it != m_stack.end()) {
      if (counter % 10000 == 0) std::cout << counter / 10000 << " " << nActive << "\n";
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

  void dumpStack(std::string filename) {
    utils::OutputFile out(filename.c_str());
    for (const auto& particle : m_stack) {
      if (particle.getPid() == proton) out << particle << "\n";
    }
  }

  double getObservedEnergy() {
    assert(m_stack.size() == 1 && m_stack[0].getRedshift() < 1e-20);
    return m_stack[0].getGamma() * SI::protonMassC2;
  }

  virtual ~Evolutor() = default;
};

void testSingleProtonEvolution() {
  RandomNumberGenerator rng = utils::RNG<double>(Seed(96));
  const double zMax = 1.;
  utils::OutputFile out("test_traj_z1.0.txt");
  for (double E = 1e17 * SI::eV; E < 1e23 * SI::eV; E *= 1.2) {
    Evolutor evolutor(rng);
    evolutor.buildParticleStack(Redshift(zMax), LorentzFactor(E / SI::protonMassC2), 1);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_evolution_1_1e17.txt");
    out << std::scientific << E / SI::eV << " " << evolutor.getObservedEnergy() / SI::eV << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    // testSingleProtonEvolution();
    {
      RandomNumberGenerator rng = utils::RNG<double>(Seed(96));
      utils::Timer timer("timer for first test");
      Evolutor evolutor(rng);
      evolutor.buildCosmologicalParticleStack(2.7, 0., 1.0, 1e3);
      evolutor.buildPhotonFields();
      evolutor.buildContinuousLosses();
      evolutor.buildStochasticInteractions();
      evolutor.run("test_proton_cosmology.txt");
      evolutor.dumpStack("test_new_spectrum_z1.0_m0.txt");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
