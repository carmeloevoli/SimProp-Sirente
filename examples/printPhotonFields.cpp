#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    const auto T_CMB = cosmo::Planck2018().tCmb;
    const auto cmb = photonfield::CMB(T_CMB);
    const auto ebl = photonfield::Dominguez2011PhotonField();
    const auto ePhoton = utils::LogAxis(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
    utils::OutputFile out("test_photonfields.txt");
    const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
    out << "# energy [eV] - lambda [mu] - \n";
    out << std::scientific;
    for (auto E : ePhoton) {
      out << E / SI::eV << "\t";
      out << energyToWavelenght(E) / SI::micron << "\t";
      out << pow2(E) * cmb.density(E) / units << "\t";
      out << pow2(E) * ebl.density(E) / units << "\t";
      out << pow2(E) * cmb.density(E, 1.) / pow3(2.) / units << "\t";
      out << pow2(E) * ebl.density(E, 1.) / pow3(2.) / units << "\t";
      out << pow2(E) * cmb.density(E, 3.) / pow3(4.) / units << "\t";
      out << pow2(E) * ebl.density(E, 3.) / pow3(4.) / units << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
