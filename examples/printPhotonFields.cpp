#include "simprop.h"

using namespace simprop;

std::pair<double, double> getDensities(const photonfields::CMB& cmb,
                                       const photonfields::PhotonField& ebl, double E, double z) {
  return {cmb.density(E, z), ebl.density(E, z)};
}

int main() {
  try {
    utils::startup_information();
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
      {
        auto n = getDensities(cmb, ebl, E, 0.);
        out << pow2(E) * n.first / units << "\t";
        out << pow2(E) * n.second / units << "\t";
      }
      {
        auto n = getDensities(cmb, ebl, E, 1.);
        out << pow2(E) * n.first / units << "\t";
        out << pow2(E) * n.second / units << "\t";
      }
      {
        auto n = getDensities(cmb, ebl, E, 2.);
        out << pow2(E) * n.first / units << "\t";
        out << pow2(E) * n.second / units << "\t";
      }
      {
        auto n = getDensities(cmb, ebl, E, 3.);
        out << pow2(E) * n.first / units << "\t";
        out << pow2(E) * n.second / units << "\t";
      }
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
