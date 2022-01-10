#include "simprop.h"

using namespace simprop;

int main() {
  try {
    // utils::startup_information();
    // auto cosmology = cosmo::Planck2018();
    // auto redshifts = utils::LinAxis(0, 20, 200);
    // utils::OutputFile out("test_cosmology.txt");
    // out << "# redshift - Hubble Time [Gyr]\n";
    // out << std::scientific;
    // for (auto z : redshifts) {
    //   out << z << "\t";
    //   out << cosmology.hubbleTime(z) / SI::Gyr << "\n";
    // }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
