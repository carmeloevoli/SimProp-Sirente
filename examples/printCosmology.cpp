#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("timer for Cosmology");
    auto cosmology = cosmo::Planck2018();
    LOGD << "Universe Age : " << cosmology.redshift2UniverseAge(0.) / SI::Gyr;
    auto redshifts = utils::LinAxis(0, 10, 10000);
    utils::OutputFile out("test_cosmology.txt");
    out << "# redshift - Hubble Time [Gyr]\n";
    out << std::scientific;
    for (auto z : redshifts) {
      out << z << "\t";
      out << 1. / cosmology.hubbleRate(z) / SI::Gyr << "\t";
      out << cosmology.redshift2ComovingDistance(z) / SI::Mpc << "\t";
      out << cosmology.redshift2LuminosityDistance(z) / SI::Mpc << "\t";
      out << cosmology.redshift2LightTravelDistance(z) / SI::Mpc << "\t";
      out << cosmology.redshift2UniverseAge(z) / SI::Gyr << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
