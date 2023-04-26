#include "simprop.h"

using namespace simprop;

void makeEnergyLossesTables() {
  solutions::Beniamino b;
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);

  auto E = utils::LogAxis(1e17 * SI::eV, 1e24 * SI::eV, 71);
  utils::OutputFile out("SimProp_proton_losses.txt");
  out << "#log10(E) [eV] - log10(beta_0) [1/yr] - log10(db_0/dE) [1/yr] - log10(beta_a) [1/yr]\n";
  for (const auto &E_i : E) {
    out << std::scientific << std::log10(E_i / SI::eV) << "\t";
    // out << std::log10(b.beta(E_i) * SI::year) << "\t";
    // out << std::log10(b.dbdE(E_i) * SI::year) << "\t";
    // out << std::log10(adiabatic.beta(proton, E_i / SI::protonMassC2) * SI::year) << "\t";
    out << "\n";
  }
}

void testCharacteristics() {
  solutions::Beniamino b;
  {
    utils::OutputFile out("SimProp_characteristics_vs_redshift.txt");
    auto z = utils::LogAxis(1e-4, 5., 1000);
    out << "# z - 10^17 eV - 10^18 eV - 10^19 eV - 10^20 eV - 10^21 eV\n";
    for (const auto &z_i : z) {
      out << std::scientific << z_i << "\t";
      out << b.generationEnergy(1e17 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e18 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e19 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e20 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e21 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("SimProp_characteristics_vs_energy.txt");
    auto E = utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
    out << "# E - 0.05 - 0.5 - 1 - 2 - 3 - 5\n";
    for (const auto &E_i : E) {
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 0.05, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 0.5, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 1., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 2., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 3., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 5., 1e-5) / SI::eV << "\t";
      out << "\n";
    }
  }
}

void testJacobian() {
  solutions::Beniamino b;
  {
    utils::OutputFile out("SimProp_jacobian_vs_redshift.txt");
    auto z = utils::LogAxis(1e-4, 5., 1000);
    out << "# z\n";
    for (const auto &z_i : z) {
      out << std::scientific << z_i << "\t";
      out << b.dilationFactor(1e17 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e18 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e19 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e20 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e21 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("SimProp_jacobian_vs_energy.txt");
    auto E = utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
    out << "# E - 0.05 - 0.5 - 1 - 2 - 3 - 5\n";
    for (const auto &E_i : E) {
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.dilationFactor(E_i, 0., 0.05, 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 0.5, 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 1., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 2., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 3., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 5., 1e-6) << "\t";
      out << "\n";
    }
  }
}

void testSpectrum() {
  solutions::Beniamino b;
  auto E = utils::LogAxis(1e17 * SI::eV, 1e23 * SI::eV, 8 * 4);
  {
    utils::OutputFile out("SimProp_proton_spectrum.txt");
    const double zMax = 6.;
    const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i, 0., zMax) / units << "\t";
      out << b.computeFlux(E_i, 1., zMax) / units << "\t";
      out << b.computeFlux(E_i, 2., zMax) / units << "\t";
      out << b.computeFlux(E_i, 3., zMax) / units << "\t";
      out << b.computeFlux(E_i, 4., zMax) / units << "\t";
      out << b.computeFlux(E_i, 5., zMax) / units << "\t";
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for analytical solution");
    // makeEnergyLossesTables();
    // testCharacteristics();
    // testJacobian();
    testSpectrum();
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
