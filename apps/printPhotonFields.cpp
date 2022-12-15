#include "simprop.h"

using namespace simprop;

void plot_local_photonfields() {
  std::vector<std::shared_ptr<photonfields::PhotonField> > phFields;
  phFields.push_back(std::make_shared<photonfields::CMB>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011UpperPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021UpperPhotonField>());
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("test_local_photonfields.txt");
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  out << "# energy [eV] - lambda [mu] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto phField : phFields) out << pow2(E) * phField->density(E) / units << "\t";
    out << "\n";
  }
}

void plot_photon_integral() {
  const auto cmb = photonfields::CMB();
  const auto ebl = photonfields::Dominguez2011PhotonField();
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("test_photon_integral.txt");
  const auto units = 1. / pow2(SI::eV) / SI::m3;
  out << "# energy [eV] - I_gamma [] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    out << E * cmb.I_gamma(E) / units << "\t";
    out << E * ebl.I_gamma(E) / units << "\t";
    out << "\n";
  }
}

void plot_evolution_photonfields() {
  const auto cmb = photonfields::CMB();
  LOGD << "CMB photon range : " << cmb.getMinPhotonEnergy() / SI::eV << " "
       << cmb.getMaxPhotonEnergy() / SI::eV;
  const auto ebl = photonfields::Dominguez2011PhotonField();
  LOGD << "EBL photon range : " << ebl.getMinPhotonEnergy() / SI::eV << " "
       << ebl.getMaxPhotonEnergy() / SI::eV;
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("test_photonfields.txt");
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  out << "# energy [eV] - lambda [mu] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    out << pow2(E) * cmb.density(E) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.1) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.2) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.5) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.75) / units << "\t";
    out << pow2(E) * ebl.density(E, 1.0) / units << "\t";
    out << pow2(E) * ebl.density(E, 1.5) / units << "\t";
    out << pow2(E) * ebl.density(E, 2.0) / units << "\t";
    out << pow2(E) * ebl.density(E, 2.5) / units << "\t";
    out << pow2(E) * ebl.density(E, 3.0) / units << "\t";
    out << "\n";
  }
}

void plot_tau() {
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  const auto cosmology = std::make_shared<cosmo::Cosmology>(0.7, 0., 0.3 * pow2(0.7), 0.7);

  const auto tau = core::OpticalDepth(cosmology, ebl);
  const auto eGamma = utils::LogAxis<double>(1e-2 * SI::TeV, 1e2 * SI::TeV, 4 * 16);
  utils::OutputFile out("test_optical_depth.txt");
  out << "# energy [eV] - tau\n";
  out << std::scientific;
  for (auto E : eGamma) {
    std::cout << E / SI::TeV << "\n";
    out << E / SI::TeV << "\t";
    out << tau.get(E, 0.1) << "\t";
    out << tau.get(E, 0.3) << "\t";
    out << tau.get(E, 0.6) << "\t";
    out << tau.get(E, 1.0) << "\t";
    out << tau.get(E, 2.0) << "\t";
    out << "\n";
  }
}

void map_tau() {
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  const auto cosmology = std::make_shared<cosmo::Cosmology>(0.7, 0., 0.3 * pow2(0.7), 0.7);
  const auto tau = core::OpticalDepth(cosmology, ebl);
  const auto logEnergyGamma = utils::LinAxis<double>(-2., 2., 101);
  const auto redshifts = utils::LinAxis<double>(0, 2, 101);
  utils::OutputFile out("map_optical_depth_cmb.txt");
  out << "# logEnergy - z - tau\n";
  out << std::scientific;
  for (auto logE : logEnergyGamma) {
    auto E = std::pow(10., logE) * SI::TeV;
    for (auto z : redshifts) {
      std::cout << logE << "\t" << z << "\n";
      out << logE << "\t" << z << "\t" << tau.get(E, z) << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    plot_local_photonfields();
    // plot_photon_integral();
    // plot_evolution_photonfields();
    // map_tau();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
