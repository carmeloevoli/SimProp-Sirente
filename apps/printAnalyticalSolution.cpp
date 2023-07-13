#include "simprop.h"

using namespace simprop;

void printLosses() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);

  auto cmb = std::make_shared<photonfields::CMB>();
  auto pair_cmb = losses::PairProductionLosses(cmb);
  auto phpion_cmb = losses::PhotoPionContinuousLosses(cmb);

  const auto gammaAxis = utils::LogAxis<double>(1e7, 1e13, 6 * 32);

  utils::OutputFile out("SimProp_proton_lambda.txt");
  out << std::scientific;
  for (auto Gamma : gammaAxis) {
    out << Gamma * SI::protonMassC2 / SI::eV << "\t";
    out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair_cmb.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_cmb.beta(proton, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }
}

void makeEnergyLossesTables(const solutions::Beniamino &b) {
  auto E = utils::LogAxis(1e15 * SI::eV, 1e24 * SI::eV, 71);
  utils::OutputFile out("SimProp_proton_losses.txt");
  out << "#log10(E) [eV] - log10(beta_0) [1/yr] - log10(db_0/dE) [1/yr]\n";
  for (const auto &E_i : E) {
    out << std::scientific << std::log10(E_i / SI::eV) << "\t";
    out << std::log10(b.beta(E_i) * SI::year) << "\t";
    out << std::log10(b.dbdE(E_i) * SI::year) << "\t";
    out << "\n";
  }
}

void testCharacteristics(const solutions::Beniamino &b) {
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

void testJacobian(const solutions::Beniamino &b) {
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

void testSpectrum(const solutions::Beniamino &b, const std::string &filename) {
  auto E = utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 16 * 4);
  {
    utils::OutputFile out(filename);
    const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.computeFlux(E_i, 0.) / units << "\t";
      out << b.computeFlux(E_i, 1.) / units << "\t";
      out << b.computeFlux(E_i, 2.) / units << "\t";
      out << b.computeFlux(E_i, 3.) / units << "\t";
      out << b.computeFlux(E_i, 4.) / units << "\t";
      out << b.computeFlux(E_i, 5.) / units << "\t";
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for analytical solution");
    printLosses();

    auto cosmology = std::make_shared<cosmo::Planck2018>();
    auto cmb = std::make_shared<photonfields::CMB>();
    std::vector<std::shared_ptr<losses::ContinuousLosses>> losses{
        std::make_shared<losses::PairProductionLosses>(cmb),
        std::make_shared<losses::PhotoPionContinuousLosses>(cmb)};

    {
      auto b = solutions::Beniamino({2.01, 0, -1, 10.}, cosmology, losses).doCaching();
      testSpectrum(b, "SimProp_proton_spectrum_2.0_0.txt");
    }
    {
      auto b = solutions::Beniamino({2.3, -3, -1, 10.}, cosmology, losses).doCaching();
      testSpectrum(b, "SimProp_proton_spectrum_2.3_-3.txt");
    }
    {
      auto b = solutions::Beniamino({2.3, 3, -1, 10.}, cosmology, losses).doCaching();
      testSpectrum(b, "SimProp_proton_spectrum_2.3_3.txt");
    }
    {
      auto b = solutions::Beniamino({2.6, 0, -1, 10.}, cosmology, losses).doCaching();
      testSpectrum(b, "SimProp_proton_spectrum_2.6_0.txt");
    }
    {
      auto b = solutions::Beniamino({2.3, 0, -1, 10.}, cosmology, losses).doCaching();
      testSpectrum(b, "SimProp_proton_spectrum_2.3_0.txt");
    }
    // testCharacteristics(b);
    // testJacobian(b);
    // testSpectrum(b);
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
