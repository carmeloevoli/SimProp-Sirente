#include "simprop.h"

using namespace simprop;

void testSingleParticleEvolution(double zMax, std::string filename) {
  RandomNumberGenerator rng = utils::RNG<double>(69);
  auto cmb = std::make_shared<photonfields::CMB>();
  auto cosmo = std::make_shared<cosmo::Cosmology>();
  auto sim = evolutors::SingleProtonEvolutor(rng);
  sim.addCosmology(cosmo);
  sim.addLosses({std::make_shared<losses::AdiabaticContinuousLosses>(cosmo),
                 std::make_shared<losses::PairProductionLosses>(cmb),
                 std::make_shared<losses::PhotoPionContinuousLosses>(cmb)});
  utils::OutputFile out(filename);
  auto energyAxis = utils::LogAxis(1e16 * SI::eV, 1e23 * SI::eV, 7 * 4);
  size_t N = 1;
  for (auto& E : energyAxis) {
    auto builder = SingleParticleBuilder(proton, {E / SI::protonMassC2, zMax}, N);
    auto stack = builder.build();
    sim.run(stack);
    out << std::scientific << E / SI::eV << " ";
    out << stack[0].getGamma() * SI::protonMassC2 / SI::eV << " ";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    testSingleParticleEvolution(5.0, "SimProp_protontraj_z5.0.txt");
    testSingleParticleEvolution(3.0, "SimProp_protontraj_z3.0.txt");
    testSingleParticleEvolution(2.0, "SimProp_protontraj_z2.0.txt");
    testSingleParticleEvolution(1.0, "SimProp_protontraj_z1.0.txt");
    testSingleParticleEvolution(0.5, "SimProp_protontraj_z0.5.txt");
    testSingleParticleEvolution(.05, "SimProp_protontraj_z0.05.txt");
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}