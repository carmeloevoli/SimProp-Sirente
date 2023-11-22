#include "simprop.h"

using namespace simprop;

// #define VERYLARGEENERGY (1e25 * SI::eV)

// auto IsActive = [](const Particle& p) {
//   constexpr double minPropagatingGamma = 1e7;
//   return (p.isNucleus() && p.isActive() && p.getRedshift() > 1e-20 &&
//           p.getGamma() > minPropagatingGamma);
// };

// class Evolutor {
//  protected:
//   const double deltaGammaCritical = 0.1;
//   RandomNumberGenerator& m_rng;
//   ParticleStack m_stack;
//   std::shared_ptr<cosmo::Cosmology> m_cosmology;
//   std::shared_ptr<photonfields::PhotonField> m_cmb;
//   std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;

//  public:
//   Evolutor(RandomNumberGenerator& rng) : m_rng(rng) {
//     m_cosmology = std::make_shared<cosmo::Planck2018>();
//   }

//   inline const ParticleStack& getStack() const { return m_stack; }

//   void buildParticleStack(double z, double Gamma, size_t N = 1) {
//     auto builder = SingleParticleBuilder(proton, {Gamma, z}, N);
//     m_stack = builder.build();
//   }

// void run(std::string filename) {
//   utils::OutputFile out(filename.c_str());
//   // size_t nActive = std::count_if(m_stack.begin(), m_stack.end(), IsActive);
//   size_t iniSize = m_stack.size();
//   const double initEnergy = sumEnergy();
//   size_t counter = 0;
//   auto it = m_stack.begin();
//   while (it != m_stack.end()) {
//     if (counter % iniSize == 0) {
//       LOGD << counter / iniSize << "\t" << sumEnergy() / initEnergy << "\t" << m_stack.size();
//     }
//     const auto nowRedshift = it->getRedshift();
//     const auto Gamma = it->getGamma();
//     const auto dz_c = computeLossesRedshiftInterval(*it);
//     assert(dz_c > 0. && dz_c <= nowRedshift);
//     const auto deltaGamma = computeDeltaGamma(*it, dz_c);
//     it->getNow() = {nowRedshift - dz_c, Gamma * (1. - deltaGamma)};
//     it = std::find_if(it, m_stack.end(), IsActive);
//     counter++;
//   }
// }  // run()

// double getObservedEnergy() {  // TODO remove this
//   assert(m_stack.size() == 1 && m_stack[0].getRedshift() < 1e-20);
//   return m_stack[0].getGamma() * SI::protonMassC2;
// }

//   double sumEnergy() {
//     double value = 0;
//     std::for_each(m_stack.begin(), m_stack.end(), [&](const auto& particle) {
//       value += (IsActive(particle)) ? particle.getEnergy() : 0.;
//     });
//     return value;
//   }

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

// void evolvePopulation() {
//   RandomNumberGenerator rng = utils::RNG<double>(69);
//   {
//     Evolutor evolutor(rng);
//     evolutor.buildCosmologicalParticleStack(2.6, 0., 1.0, 1e5);
//     evolutor.buildPhotonFields();
//     evolutor.buildContinuousLosses(true);
//     evolutor.run("test_proton_cosmology.txt");
//     evolutor.dump("test_spectrum_z1.0_m0_N1e6.txt");
//   }
// }

void testSpectrumEvolution(double zMax, std::string filename, size_t N = 100) {
  RandomNumberGenerator rng = utils::RNG<double>(69);
  auto cmb = std::make_shared<photonfields::CMB>();
  auto cosmo = std::make_shared<cosmo::Cosmology>();
  auto sim = evolutors::SingleProtonEvolutor(rng);
  sim.addCosmology(cosmo);
  sim.addLosses({std::make_shared<losses::AdiabaticContinuousLosses>(cosmo),
                 std::make_shared<losses::PairProductionLosses>(cmb)});
  sim.addInteractions({std::make_shared<interactions::PhotoPionProduction>(cmb)});

  const auto minEnergy = 1e17 * SI::eV;
  const auto maxEnergy = 1e23 * SI::eV;
  const auto slope = 2.6;
  const auto m = 0.;
  Range GammaRange = {minEnergy / SI::protonMassC2, maxEnergy / SI::protonMassC2};
  Range zRange = {0., zMax};
  auto builder = SourceEvolutionBuilder(proton, {GammaRange, zRange, slope, m}, cosmo, N);
  auto stack = builder.build(rng);
  sim.run(stack);

  utils::OutputFile out(filename);
  for (const auto& particle : stack) {
    if (particle.getPid() == proton && particle.getRedshift() < 1e-20) out << particle << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    testSpectrumEvolution(3.0, "SimProp_spectrum_z3.0.txt", 100000);
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}