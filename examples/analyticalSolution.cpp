#include "simprop.h"

using namespace simprop;

void testCharacteristics() {
  solutions::Beniamino b;
  {
    utils::OutputFile out("proton_characteristics_redshift.txt");
    auto z = utils::LogAxis(1e-4, 6., 1000);
    for (const auto &z_i : z) {
      std::cout << z_i << "\n";
      out << std::scientific << z_i << "\t";
      out << b.generationEnergy(1e17 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e18 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e19 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e20 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e21 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(1e22 * SI::eV, z_i, 1e-6) / SI::eV << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("proton_characteristics_energy.txt");
    auto E = utils::LogAxis(1e16 * SI::eV, 1e23 * SI::eV, 1000);
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0.05, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0.5, 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 1., 1e-6) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 3., 1e-6) / SI::eV << "\t";
      out << "\n";
    }
  }
}

void testJacobian() {
  solutions::Beniamino b;
  b.disablePhotoPion();
  utils::OutputFile out("proton_dilationFactor.txt");
  auto z = utils::LogAxis(1e-4, 1., 100);
  for (auto it = z.begin(); it != z.end(); it++) {
    auto z_i = *it;
    std::cout << z_i << "\n";
    out << std::scientific << z_i << "\t";
    out << b.dilationFactor(1e17 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e19 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e21 * SI::eV, z_i) << "\t";
    out << b.dilationFactor(1e23 * SI::eV, z_i) << "\t";
    out << "\n";
  }
}

void printSpectrum() {
  solutions::Beniamino b;
  {
    b.setMaxRedshift(0.1);
    b.setSlope(2.7);
    b.setSourceEvolution(0.);
    // b.disablePhotoPion();
    utils::OutputFile out("proton_spectrum_zmax0.1_m0.txt");
    auto E = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 5 * 16);
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i) << "\t";
      out << b.computeFluxUnm(E_i) << "\t";
      out << "\n";
    }
  }
  {
    b.setMaxRedshift(0.5);
    b.setSlope(2.7);
    b.setSourceEvolution(0.);
    // b.disablePhotoPion();
    utils::OutputFile out("proton_spectrum_zmax0.5_m0.txt");
    auto E = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 5 * 16);
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i) << "\t";
      out << b.computeFluxUnm(E_i) << "\t";
      out << "\n";
    }
  }
  {
    b.setMaxRedshift(1.0);
    b.setSlope(2.7);
    b.setSourceEvolution(0.);
    // b.disablePhotoPion();
    utils::OutputFile out("proton_spectrum_zmax1.0_m0.txt");
    auto E = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 5 * 16);
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i) << "\t";
      out << b.computeFluxUnm(E_i) << "\t";
      out << "\n";
    }
  }
  {
    b.setMaxRedshift(3.0);
    b.setSlope(2.7);
    b.setSourceEvolution(0.);
    // b.disablePhotoPion();
    utils::OutputFile out("proton_spectrum_zmax3.0_m0.txt");
    auto E = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 5 * 16);
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i) << "\t";
      out << b.computeFluxUnm(E_i) << "\t";
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
    printSpectrum();
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
