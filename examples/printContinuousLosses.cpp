// Copyright 2023 SimProp-dev [MIT License]
#include "simprop.h"

using namespace simprop;

void plot_pair_on_fields() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);

  std::vector<std::shared_ptr<photonfields::PhotonField>> phFields{
      std::make_shared<photonfields::CMB>(),
      std::make_shared<photonfields::Saldana2021PhotonField>()};
  auto pair_full = losses::PairProductionLosses(phFields);

  auto cmb = std::make_shared<photonfields::CMB>();
  auto pair_cmb = losses::PairProductionLosses(cmb);

  auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  auto pair_ebl = losses::PairProductionLosses(ebl);

  auto bggLosses = losses::BGG2006ContinuousLosses();

  const auto gammaAxis = utils::LogAxis<double>(1e7, 1e13, 6 * 32);

  utils::OutputFile out("SimProp_losses_pair.txt");
  out << std::scientific;
  for (auto Gamma : gammaAxis) {
    out << Gamma << "\t";
    out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / bggLosses.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair_cmb.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair_ebl.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair_full.beta(proton, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }
}

void plot_pair_nuclei() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
  auto bggLosses = losses::BGG2006ContinuousLosses();

  std::vector<std::shared_ptr<photonfields::PhotonField>> phFields{
      std::make_shared<photonfields::CMB>(),
      std::make_shared<photonfields::Dominguez2011PhotonField>()};
  auto pair = losses::PairProductionLosses(phFields);

  const auto gammaAxis = utils::LogAxis<double>(1e7, 1e13, 6 * 32);
  utils::OutputFile out("SimProp_losses_pair_nuclei.txt");

  out << std::scientific;
  for (auto Gamma : gammaAxis) {
    out << Gamma << "\t";
    out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / bggLosses.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / bggLosses.beta(O16, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / bggLosses.beta(Fe56, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair.beta(O16, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / pair.beta(Fe56, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }
}

void plot_photopion_on_fields() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);

  std::vector<std::shared_ptr<photonfields::PhotonField>> phFields{
      std::make_shared<photonfields::CMB>(),
      std::make_shared<photonfields::Saldana2021PhotonField>()};
  auto phpion_full = losses::PhotoPionContinuousLosses(phFields);

  auto cmb = std::make_shared<photonfields::CMB>();
  auto phpion_cmb = losses::PhotoPionContinuousLosses(cmb);

  auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  auto phpion_ebl = losses::PhotoPionContinuousLosses(ebl);

  const auto gammaAxis = utils::LogAxis<double>(1e7, 1e13, 6 * 32);

  utils::OutputFile out("SimProp_losses_photopion.txt");
  out << std::scientific;
  for (auto Gamma : gammaAxis) {
    out << Gamma << "\t";
    out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_cmb.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_ebl.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_full.beta(proton, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }
}

void plot_photopion_nuclei() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
  std::vector<std::shared_ptr<photonfields::PhotonField>> phFields{
      std::make_shared<photonfields::CMB>(),
      std::make_shared<photonfields::Saldana2021PhotonField>()};
  auto phpion_full = losses::PhotoPionContinuousLosses(phFields);

  auto cmb = std::make_shared<photonfields::CMB>();
  auto phpion_cmb = losses::PhotoPionContinuousLosses(cmb);

  const auto gammaAxis = utils::LogAxis<double>(1e8, 1e13, 6 * 32);
  utils::OutputFile out("SimProp_losses_photopion_nuclei.txt");
  out << std::scientific;
  for (auto Gamma : gammaAxis) {
    out << Gamma << "\t";
    out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_cmb.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_cmb.beta(O16, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_cmb.beta(Fe56, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_full.beta(proton, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_full.beta(O16, Gamma) / SI::Mpc << "\t";
    out << SI::cLight / phpion_full.beta(Fe56, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }
}

// void plot_pair_redshift() {
//   auto cosmology = std::make_shared<cosmo::Planck2018>();
//   auto adiabatic = losses::AdiabaticContinuousLosses(cosmology);
//   auto bggLosses = losses::BGG2006ContinuousLosses();
//   auto cmb = std::make_shared<photonfields::CMB>();
//   auto pair = losses::PairProductionLosses(cmb);

//   const auto gammaAxis = utils::LogAxis<double>(1e7, 1e13, 6 * 32);
//   utils::OutputFile out("test_redshift_pair_losses.txt");

//   out << std::scientific;
//   for (auto Gamma : gammaAxis) {
//     out << Gamma << "\t";
//     out << SI::cLight / adiabatic.beta(proton, Gamma) / SI::Mpc << "\t";
//     out << SI::cLight / bggLosses.beta(proton, Gamma) / SI::Mpc << "\t";
//     out << SI::cLight / bggLosses.beta(proton, Gamma, 1.0) / SI::Mpc << "\t";
//     out << SI::cLight / bggLosses.beta(proton, Gamma, 3.0) / SI::Mpc << "\t";
//     out << SI::cLight / bggLosses.beta(proton, Gamma, 6.0) / SI::Mpc << "\t";
//     out << SI::cLight / pair.beta(proton, Gamma) / SI::Mpc << "\t";
//     out << SI::cLight / pair.beta(proton, Gamma, 1.0) / SI::Mpc << "\t";
//     out << SI::cLight / pair.beta(proton, Gamma, 3.0) / SI::Mpc << "\t";
//     out << SI::cLight / pair.beta(proton, Gamma, 6.0) / SI::Mpc << "\t";
//     out << "\n";
//   }
// }

void plot_photopion_inelasticity() {
  double eps_th = SI::pionMassC2 + SI::pionMassC2 * SI::pionMassC2 / 2. / SI::protonMassC2;
  utils::OutputFile out("SimProp_losses_photopion_inelasticity.txt");
  out << std::scientific;
  for (double eps = eps_th; eps < 10.0 * SI::GeV; eps += 0.01 * SI::GeV) {
    out << eps / SI::GeV << "\t";
    auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * eps;
    out << s / SI::GeV2 << "\t";
    out << losses::inelasticity(eps) << "\t";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for print losses");
    plot_pair_on_fields();
    plot_pair_nuclei();
    plot_photopion_on_fields();
    plot_photopion_nuclei();
    // simprop::plot_pair_redshift();
    plot_photopion_inelasticity();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
