#include "simprop.h"
#include "simprop/utils/misc.h"

using namespace simprop;

int main() {
  try {
    log::startup_information();
    const auto cmb = photonfield::CMB();
    const auto ebl = photonfield::Dominguez2011PhotonField();
    const auto ePhoton = utils::LogAxis(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
    utils::OutputFile out("test_EBL.txt");
    const auto units = SI::nW / utils::pow<2>(SI::meter) / SI::sr;
    out() << std::scientific;
    for (auto E : ePhoton) {
      out() << E / SI::eV << "\t";
      out() << energyToWavelenght(E) / SI::micron << "\t";
      const auto n_cmb = cmb.getPhotonDensity(E);
      out() << SI::cOver4pi * utils::pow<2>(E) * n_cmb / units << "\t";
      const auto n_ebl = ebl.getPhotonDensity(E);
      out() << SI::cOver4pi * utils::pow<2>(E) * n_ebl / units << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
