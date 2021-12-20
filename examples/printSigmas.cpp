#include "simprop.h"

using namespace simprop;

int main() {
  try {
    log::startup_information();
    {
      auto sigmaPpp = interactions::PhotoPionProduction();
      auto sAxis = utils::LogAxis(0.1 * SI::GeV2, 1e12 * SI::GeV2, 1000);
      utils::OutputFile out("test_ppp_sigma.txt");
      out() << std::scientific;
      for (auto s : sAxis) {
        out() << s / SI::GeV2 << "\t";
        out() << sigmaPpp.getAtS(s) / SI::mbarn << "\n";
      }
    }

    {
      auto sigmaPd = interactions::PsbPhotoDisintegration();
      auto eAxis = utils::LinAxis(0.1 * SI::MeV, 1e2 * SI::MeV, 1000);
      utils::OutputFile out("test_pd_sigma.txt");
      out() << std::scientific;
      for (auto E : eAxis) {
        out() << E / SI::MeV << "\t";
        out() << sigmaPd.getSigma(Ca40, E) / SI::mbarn << "\n";
      }
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}