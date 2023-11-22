#include "simprop.h"

using namespace simprop;

void testSpectrumEvolution(double zMax, std::string filename, size_t N = 100) {
  RandomNumberGenerator rng = utils::RNG<double>(69);
  auto cmb = std::make_shared<photonfields::CMB>();
  auto cosmo = std::make_shared<cosmo::Cosmology>();
  auto sim = evolutors::SingleProtonEvolutor(rng);
  sim.addCosmology(cosmo);
  auto adiabatic = std::make_shared<losses::AdiabaticContinuousLosses>(cosmo);
  auto pp = std::make_shared<losses::PairProductionLosses>(cmb);
  pp->doCaching();
  sim.addLosses({adiabatic, pp});
  auto ppp = std::make_shared<interactions::PhotoPionProductionSophia>(cmb);
  ppp->doCaching();
  sim.addInteractions({ppp});

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
    testSpectrumEvolution(3.0, "SimProp_spectrum_z3.0_sophia.txt", 1000000);
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}