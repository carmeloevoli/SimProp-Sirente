#include "simprop.h"

using namespace simprop;

void plot_photonfields_local(std::string filename) {
  std::vector<std::shared_ptr<photonfields::PhotonField> > phFields;
  phFields.push_back(std::make_shared<photonfields::CMB>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011UpperPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021UpperPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Gilmore2012PhotonField>());
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out(filename);
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  out << "# energy [eV] - lambda [mu] - density []\n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto phField : phFields) out << pow2(E) * phField->density(E) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_at_z(const photonfields::PhotonField& phField, std::string filename) {
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  const auto cmb = std::make_shared<photonfields::CMB>();

  utils::OutputFile out(filename);
  out << "# eps [eV] - photon density\n";
  out << std::scientific;
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  std::vector<double> redshifts{0., 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto z : redshifts)
      out << pow2(E) * phField.density(E, z) / std::pow(1. + z, 3.) / units << "\t";
    for (auto z : redshifts)
      out << pow2(E) * cmb->density(E, z) / std::pow(1. + z, 3.) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_Igamma() {
  std::vector<std::shared_ptr<photonfields::PhotonField> > phFields;
  phFields.push_back(std::make_shared<photonfields::CMB>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Gilmore2012PhotonField>());
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("test_photon_Igamma.txt");
  const auto units = 1. / SI::eV / SI::m3;
  out << "# energy [eV] - wavelength [mu] - I_gamma [eV^-1 m^-3] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto phField : phFields) out << E * phField->I_gamma(E) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_density() {
  std::vector<std::shared_ptr<photonfields::PhotonField> > phFields;
  phFields.push_back(std::make_shared<photonfields::CMB>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Dominguez2011UpperPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021PhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021LowerPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Saldana2021UpperPhotonField>());
  phFields.push_back(std::make_shared<photonfields::Gilmore2012PhotonField>());
  const auto redshifts = utils::LinAxis<double>(0, 6, 1000);
  utils::OutputFile out("test_photon_density.txt");
  out << "# z - photon density at 0.01 eV - photon density at 1 eV\n";
  out << std::scientific;
  auto eps = 1e-2 * SI::eV;
  auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  for (auto z : redshifts) {
    out << z << "\t";
    for (auto phField : phFields)
      out << eps * eps * phField->density(eps, z) / std::pow(1. + z, 3.) / units << "\t";
    for (auto phField : phFields)
      out << phField->computeIntegratedDensity(z) / std::pow(1. + z, 3.) / (1. / SI::m3) << "\t";
    out << "\n";
  }
}

void plot_tau() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  const auto cosmology = std::make_shared<cosmo::Cosmology>(0.7, 0., 0.3 * pow2(0.7), 0.7);
  const auto tau_cmb = core::OpticalDepth(cosmology, cmb);
  const auto tau_ebl = core::OpticalDepth(cosmology, ebl);
  const auto eGamma = utils::LogAxis<double>(1e-2 * SI::TeV, 1e2 * SI::TeV, 4 * 8);
  const auto z = std::vector<double>({0.1, 0.5, 1.0, 3.0, 5.0});
  utils::OutputFile out("test_optical_depth.txt");
  out << "# energy [TeV] - tau [z=0.1,0.5,1,2,3]\n";
  out << std::scientific;
  for (auto E : eGamma) {
    std::cout << E / SI::TeV << "\n";
    out << E / SI::TeV << "\t";
    for (auto z_i : z) out << tau_cmb.get(E, z_i) << "\t";
    for (auto z_i : z) out << tau_ebl.get(E, z_i) << "\t";
    out << "\n";
  }
}

void map_tau() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  const auto cosmology = std::make_shared<cosmo::Cosmology>(0.7, 0., 0.3 * pow2(0.7), 0.7);
  const auto tau = core::OpticalDepth(cosmology, {cmb, ebl});
  const auto logEnergyGamma = utils::LinAxis<double>(13., 22., 101);
  const auto logDistances = utils::LinAxis<double>(-3., 2., 101);
  utils::OutputFile out("map_optical_depth_ebl.txt");
  out << "# logEnergy - logDistance - tau\n";
  out << std::scientific;
  for (auto logE : logEnergyGamma) {
    auto E = std::pow(10., logE) * SI::eV;
    for (auto logd : logDistances) {
      auto d = std::pow(10., logd) * SI::Mpc;
      auto z = cosmology->comovingDistance2Redshift(d);
      std::cout << logE << "\t" << logd << "\t" << z << "\n";
      out << logE << "\t" << logd << "\t" << tau.get(E, z) << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    // plot_photonfields_local("test_local_photonfields.txt");
    // plot_photonfields_at_z(photonfields::Saldana2021PhotonField(),
    //                        "test_photonfields_Saldana2021.txt");
    // plot_photonfields_at_z(photonfields::Saldana2021LowerPhotonField(),
    //                        "test_photonfields_Saldana2021_lower.txt");
    // plot_photonfields_at_z(photonfields::Saldana2021UpperPhotonField(),
    //                        "test_photonfields_Saldana2021_upper.txt");
    // plot_photonfields_at_z(photonfields::Dominguez2011PhotonField(),
    //                        "test_photonfields_Dominguez2011.txt");
    // plot_photonfields_at_z(photonfields::Dominguez2011LowerPhotonField(),
    //                        "test_photonfields_Dominguez2011_lower.txt");
    // plot_photonfields_at_z(photonfields::Dominguez2011UpperPhotonField(),
    //                        "test_photonfields_Dominguez2011_upper.txt");
    // plot_photonfields_at_z(photonfields::Gilmore2012PhotonField(),
    //                        "test_photonfields_Gilmore2012.txt");
    // plot_photonfields_Igamma();
    // plot_photonfields_density();
    plot_tau();
    // map_tau();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
