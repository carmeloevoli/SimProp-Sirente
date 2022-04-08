#include "simprop.h"

using namespace simprop;

void plot_local_photonfields() {
  const auto cmb = photonfields::CMB();
  const auto ebl = photonfields::Dominguez2011PhotonField();
  const auto ebl_lo = photonfields::Dominguez2011PhotonField(photonfields::EblModel::LOWER);
  const auto ebl_up = photonfields::Dominguez2011PhotonField(photonfields::EblModel::UPPER);
  const auto ePhoton = utils::LogAxis(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("test_local_photonfields.txt");
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  out << "# energy [eV] - lambda [mu] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    out << pow2(E) * cmb.density(E) / units << "\t";
    out << pow2(E) * ebl.density(E, 0.) / units << "\t";
    out << pow2(E) * ebl_lo.density(E, 0.) / units << "\t";
    out << pow2(E) * ebl_up.density(E, 0.) / units << "\t";
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
  const auto ePhoton = utils::LogAxis(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
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

int main() {
  try {
    utils::startup_information();
    plot_local_photonfields();
    plot_evolution_photonfields();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
