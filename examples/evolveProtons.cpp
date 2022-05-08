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

  void buildCosmologicalParticleStack(size_t N = 1) {
    // parameters are taken from R. Alves Batista, et al., JCAP 2015, arXiv:1508.01824
    double maxEnergy = std::pow(10., 23) * SI::eV;
    double minEnergy = std::pow(10., 17) * SI::eV;
    double slope = 2.7;
    double evolutionIndex = 0;
    Range GammaRange = {minEnergy / SI::protonMassC2, maxEnergy / SI::protonMassC2};
    Range zRange = {0., 1.0};
    auto builder = SourceEvolutionBuilder(proton, {GammaRange, zRange, slope}, N);
    m_stack = builder.build(m_rng);
    dumpStack("initial_state.txt");
  }

  void buildPhotonFields() {
    m_cmb = std::make_shared<photonfields::CMB>();
    // m_ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  }

  void buildContinuousLosses() {
    // std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{m_cmb, m_ebl};
    // std::vector<std::shared_ptr<photonfields::PhotonField> > phFields{m_cmb};
    // m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
    //     std::make_shared<losses::PairProductionLosses>(phFields),
    //     std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology)};
    m_continuousLosses = std::vector<std::shared_ptr<losses::ContinuousLosses> >{
        std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology),
        std::make_shared<losses::BGG2002ContinuousLosses>()};
  }

  void buildStochasticInteractions() {
    auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
    m_pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, m_cmb);
  }

  // double computeStochasticRedshiftInterval(const Particle& particle, RandomNumber r) {
  //   const auto zNow = particle.getRedshift();
  //   const auto Gamma = particle.getGamma();
  //   const auto pid = particle.getPid();
  //   const auto dtdz = m_cosmology->dtdz(zNow);
  //   const auto lambda_s = std::fabs(1. / m_pppcmb->rate(pid, Gamma, zNow)) / dtdz;
  //   // TODO why to put the fabs?
  //   return -lambda_s * std::log(1. - r.get());
  // }

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
    size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
    // auto it = m_stack.begin();
    size_t counter = 0;
    // while (it != m_stack.end())
    while (nActive > 0) {
      if (counter % 1000 == 0) std::cout << counter / 1000 << "\n";
      const auto it = std::find_if(m_stack.begin(), m_stack.end(), IsActive);
      out << *it << " " << 0 << "\n";

      const auto nowRedshift = it->getRedshift();

      // auto r = m_rng();
      //  const auto dz_s = computeStochasticRedshiftInterval(*it, RandomNumber(r));
      //  assert(dz_s > 0.);

      const auto dz_c = computeLossesRedshiftInterval(*it);
      assert(dz_c > 0. && dz_c <= nowRedshift);

      // if (dz_s > dz_c || dz_s > nowRedshift) {
      const auto Gamma = it->getGamma();
      const auto dz = dz_c;
      const auto deltaGamma = computeDeltaGamma(*it, dz);
      it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
      // out << *it << " " << 0 << "\n";
      //   } else {
      //     // out << *it << " " << 1 << "\n";
      //     const auto dz = dz_s;
      //     auto finalState = m_pppcmb->finalState(*it, nowRedshift - dz, m_rng);
      //     m_stack.erase(it);
      //     m_stack.insert(m_stack.begin(), finalState.begin(), finalState.end());
      //   }

      nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
      counter++;
    }
  }  // run()

  void runNew(std::string filename) {
    utils::OutputFile out(filename.c_str());
    size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
    size_t counter = 0;
    // while (it != m_stack.end())

    const double dz = 1e-4;

    auto it = m_stack.begin();
    while (it != m_stack.end()) {
      if (counter % 1000 == 0) std::cout << counter / 1000 << " " << nActive << "\n";
      // const auto it = std::find_if(m_stack.begin(), m_stack.end(), IsActive);
      // out << *it << " " << 0 << "\n";
      const auto nowRedshift = it->getRedshift();
      const auto Gamma = it->getGamma();
      const auto deltaGamma = computeDeltaGamma(*it, dz);
      it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
      // nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
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
  const double zMax = 0.05;
  utils::OutputFile out("test_traj_z0.05.txt");
  for (double E = 1e17 * SI::eV; E < 1e23 * SI::eV; E *= 1.2) {
    utils::Timer timer("timer for Gamma = 1e10");
    Evolutor evolutor(rng);
    evolutor.buildParticleStack(Redshift(zMax), LorentzFactor(E / SI::protonMassC2), 1);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.buildStochasticInteractions();
    evolutor.runNew("test_proton_evolution_1_1e17.txt");
    out << std::scientific << E / SI::eV << " " << evolutor.getObservedEnergy() / SI::eV << "\n";
  }

  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(2.0017476695909485e+17 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e17.txt");
  // }
  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(8.707673289496823e+19 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e18.txt");
  // }
  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(8.308907186554102e+20 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e19.txt");
  // }
  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(1.837157436044067e+21 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e20.txt");
  // }
  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(5.280634760793176e+21 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e21.txt");
  // }
  // {
  //   utils::Timer timer("timer for Gamma = 1e10");
  //   Evolutor evolutor(rng);
  //   evolutor.buildParticleStack(
  //       Redshift(1.), LorentzFactor(2.704500955346506e+22 * SI::eV / SI::protonMassC2), 1);
  //   evolutor.buildPhotonFields();
  //   evolutor.buildContinuousLosses();
  //   evolutor.buildStochasticInteractions();
  //   evolutor.runNew("test_proton_evolution_1_1e22.txt");
  // }
}

int main() {
  try {
    utils::startup_information();
    // testSingleProtonEvolution();
    //   {
    //     RandomNumberGenerator rng = utils::RNG<double>(Seed(88));
    //     utils::Timer timer("timer for Gamma = 1e11");
    //     Evolutor evolutor(rng);
    //     evolutor.buildParticleStack(Redshift(1.), LorentzFactor(1e11), 1);
    //     evolutor.buildPhotonFields();
    //     evolutor.buildContinuousLosses();
    //     evolutor.buildStochasticInteractions();
    //     evolutor.run("test_proton_evolution_1_1e11_1.txt");
    //   }
    //   {
    //     RandomNumberGenerator rng = utils::RNG<double>(Seed(88));
    //     utils::Timer timer("timer for Gamma = 1e12");
    //     Evolutor evolutor(rng);
    //     evolutor.buildParticleStack(Redshift(1.), LorentzFactor(1e12), 1);
    //     evolutor.buildPhotonFields();
    //     evolutor.buildContinuousLosses();
    //     evolutor.buildStochasticInteractions();
    //     evolutor.run("test_proton_evolution_1_1e12_1.txt");
    //   }

    {
      RandomNumberGenerator rng = utils::RNG<double>(Seed(96));
      utils::Timer timer("timer for first test");
      Evolutor evolutor(rng);
      evolutor.buildCosmologicalParticleStack(1e6);
      evolutor.buildPhotonFields();
      evolutor.buildContinuousLosses();
      // evolutor.buildStochasticInteractions();
      evolutor.runNew("test_proton_cosmology.txt");
      evolutor.dumpStack("test_spectrum_z1.0.txt");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
