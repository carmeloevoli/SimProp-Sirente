#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    const auto Gammas = utils::LogAxis(1e8, 1e16, 8 * 32);
    const auto units = 1. / SI::Mpc;
    const auto cmb = std::make_shared<photonfields::CMB>();
    const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
    const auto sigma = std::make_shared<xsecs::PhotoPionProduction>();
    const auto pppcmb = interactions::pppEbl(sigma, cmb);
    const auto pppebl = interactions::pppEbl(sigma, ebl);
    utils::OutputFile out("test_rates.txt");
    out << std::scientific;
    for (auto Gamma : Gammas) {
      out << Gamma << "\t";
      out << pppcmb.rate(proton, Gamma) / SI::cLight / units << "\t";
      out << pppebl.rate(proton, Gamma) / SI::cLight / units << "\t";
      out << pppebl.rate(proton, Gamma, 1.0) / SI::cLight / units << "\t";
      out << pppebl.rate(proton, Gamma, 2.0) / SI::cLight / units << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
