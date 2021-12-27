#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto sigmaPpp = xsecs::PhotoPionProduction();
    auto eAxis = utils::LogAxis(1e7 * SI::eV, 1e21 * SI::eV, 1000);
    utils::OutputFile out("test_xsecs.txt");
    out << "#\n";
    out << std::scientific;
    for (auto eps : eAxis) {
      out << eps / SI::eV << "\t";
      out << sigmaPpp.get(proton, eps) / SI::mbarn << "\n";
    }

    // {
    //   auto sigmaPd = interactions::PsbPhotoDisintegration();
    //   auto eAxis = utils::LinAxis(0.1 * SI::MeV, 1e2 * SI::MeV, 1000);
    //   utils::OutputFile out("test_pd_sigma.txt");
    //   out() << std::scientific;
    //   for (auto E : eAxis) {
    //     out() << E / SI::MeV << "\t";
    //     out() << sigmaPd.getSigma(Ca40, E) / SI::mbarn << "\n";
    //   }
    // }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}