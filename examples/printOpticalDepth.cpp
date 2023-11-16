#include "simprop.h"

using namespace simprop;

void plot_single_tau(std::shared_ptr<photonfields::PhotonField> photonField, std::string filename) {
  const auto cosmology = std::make_shared<cosmo::Planck2018>();
  const auto tau = core::OpticalDepth(cosmology, photonField);
  const auto eGamma = utils::LogAxis<double>(1e-2 * SI::TeV, 1e5 * SI::TeV, 100);
  const auto z = std::vector<double>({0.1, 0.5, 1.0, 3.0, 5.0});
  utils::OutputFile out(filename);
  out << "# energy [TeV] - tau [z=0.1,0.5,1,2,3]\n";
  out << std::scientific;
  for (auto E : eGamma) {
    std::cout << std::scientific << E / SI::TeV << "\n";
    out << E / SI::TeV << "\t";
    for (auto z_i : z) out << tau.get(E, z_i) << "\t";
    out << "\n";
  }
}

void plot_tau() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  plot_single_tau(cmb, "SimProp_opticaldepth_cmb.txt");
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  plot_single_tau(ebl, "SimProp_opticaldepth_ebl.txt");
  const auto radio = std::make_shared<photonfields::Nitu2021RadioPhotonField>();
  plot_single_tau(radio, "SimProp_opticaldepth_radio.txt");
}

void plot_single_interaction_lenght(std::shared_ptr<photonfields::PhotonField> photonField,
                                    std::string filename) {
  const auto cosmology = std::make_shared<cosmo::Planck2018>();
  const auto tau = core::OpticalDepth(cosmology, photonField);
  const auto eGamma = utils::LogAxis<double>(SI::TeV, 1e14 * SI::TeV, 1000);
  utils::OutputFile out(filename);
  out << "# energy [TeV] - length [Mpc]\n";
  out << std::scientific;
  for (auto E : eGamma) {
    out << E / SI::TeV << "\t";
    out << tau.interactionLength(E) / SI::Mpc << "\t";
    out << "\n";
  }
}

void plot_interaction_lenght() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  plot_single_interaction_lenght(cmb, "SimProp_gammalenght_cmb.txt");
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  plot_single_interaction_lenght(ebl, "SimProp_gammalenght_ebl.txt");
  const auto radio = std::make_shared<photonfields::Nitu2021RadioPhotonField>();
  plot_single_interaction_lenght(radio, "SimProp_gammalenght_radio.txt");
}

void map_tau() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  const auto radio = std::make_shared<photonfields::Nitu2021RadioPhotonField>();
  const auto cosmology = std::make_shared<cosmo::Planck2018>();

  const auto tau = core::OpticalDepth(cosmology, {cmb, ebl, radio});
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
    // plot_tau();
    //  plot_interaction_lenght();
    map_tau();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
