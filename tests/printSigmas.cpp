#include "simprop.h"
#include "simprop/utils/misc.h"

using namespace simprop;

int main() {
  try {
    log::startup_information();
    auto photoPion = interactions::PhotoPion();
    auto sAxis = utils::LogAxis(0.1 * SI::GeV2, 1e12 * SI::GeV2, 1000);
    utils::OutputFile out("test_sigmas.txt");
    out() << std::scientific;
    for (auto s : sAxis) {
      out() << s / SI::GeV2 << "\t";
      out() << photoPion.sigma_of_s(s) / SI::mbarn << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}