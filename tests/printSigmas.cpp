#include "simprop.h"

using namespace simprop;

int main() {
  try {
    log::startup_information();
    auto sigmaPpp = interactions::PhotoPionProduction();
    auto sigmaPd = interactions::PsbPhotoDisintegration();
    auto sAxis = utils::LogAxis(0.1 * SI::GeV2, 1e12 * SI::GeV2, 1000);
    utils::OutputFile out("test_sigmas.txt");
    out() << std::scientific;
    for (auto s : sAxis) {
      out() << s / SI::GeV2 << "\t";
      out() << sigmaPpp.getAtS(s) / SI::mbarn << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}