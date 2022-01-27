#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto sigmaPpp = xsecs::PhotoPionProduction();
    {
      auto sAxis = utils::LogAxis(0.1 * pow2(SI::GeV), 1e11 * pow2(SI::GeV), 1000);
      utils::OutputFile out("test_xsecs.txt");
      out << "#\n";
      out << std::scientific;
      for (auto s : sAxis) {
        out << s / pow2(SI::GeV) << "\t";
        out << sigmaPpp.getAtS(proton, s) / SI::mbarn << "\n";
      }
    }
    {
      auto eAxis = utils::LogAxis(0.1 * SI::GeV, 1e11 * SI::GeV, 1000);
      utils::OutputFile out("test_xsecs_epsprime.txt");
      out << "#\n";
      out << std::scientific;
      for (auto e : eAxis) {
        out << e / SI::GeV << "\t";
        out << sigmaPpp.getAtEpsPrime(proton, e) / SI::mbarn << "\n";
      }
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