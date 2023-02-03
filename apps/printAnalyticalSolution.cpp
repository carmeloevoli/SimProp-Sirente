#include "simprop.h"

using namespace simprop;

void testCharacteristics() {
  solutions::Beniamino b(true);
  {
    utils::OutputFile out("proton_characteristics_redshift.txt");
    auto z = utils::LogAxis(1e-4, 5., 1000);
    out << "# z - 10^17 eV - 10^18 eV - 10^19 eV - 10^20 eV - 10^21 eV\n";
    for (const auto &z_i : z) {
      out << std::scientific << z_i << "\t";
      out << b.generationEnergy(1e17 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e18 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e19 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e20 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e21 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("proton_characteristics_energy.txt");
    auto E = utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
    out << "# E - 0.05 - 0.5 - 1 - 2 - 3\n";
    for (const auto &E_i : E) {
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0.05, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0.5, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 1., 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 2., 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 3., 1e-6) / SI::eV << "\t";
      out << "\n";
    }
  }
}

void testJacobian() {
  solutions::Beniamino b(true);
  utils::OutputFile out("proton_dilation_factor.txt");
  auto z = utils::LogAxis(1e-3, 1e1, 100);
  out << "# z\n";
  for (auto it = z.rbegin(); it != z.rend(); it++) {
    auto z_i = *it;
    out << std::scientific << z_i << "\t";
    out << b.dilationFactor(1e18 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e19 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e20 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e21 * SI::eV, z_i) << "\t";
    out << "\n";
  }
}

void testSpectrum() {
  solutions::Beniamino b(true);
  b.setSlope(2.7);
  b.setSourceEvolution(0.);
  auto E = utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 80);
  {
    utils::OutputFile out("proton_spectrum_pion_m0.txt");
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i, 0, 0.05) << "\t";
      out << b.computeFlux(E_i, 0, 0.5) << "\t";
      out << b.computeFlux(E_i, 0, 1.0) << "\t";
      out << b.computeFlux(E_i, 0, 3.0) << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("proton_unm_spectrum_pion_m0.txt");
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFluxUnm(E_i, 0.05) << "\t";
      out << b.computeFluxUnm(E_i, 0.5) << "\t";
      out << b.computeFluxUnm(E_i, 1.0) << "\t";
      out << b.computeFluxUnm(E_i, 3.0) << "\t";
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for analytical solution");
    testCharacteristics();
    testJacobian();
    // testSpectrum();
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
