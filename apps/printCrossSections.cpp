#include "simprop.h"

using namespace simprop;

void plot_pionproduction() {
  auto sigmaPpp = xsecs::PhotoPionProductionXsec();
  {
    auto sAxis = utils::LogAxis(0.1 * pow2(SI::GeV), 1e11 * pow2(SI::GeV), 10000);
    utils::OutputFile out("test_xsecs.txt");
    out << "# s [GeV^2] - sigma [mbarn] - phi []\n";
    out << std::scientific;
    for (auto s : sAxis) {
      out << s / SI::GeV2 << "\t";
      out << sigmaPpp.getAtS(s) / SI::mbarn << "\t";
      out << sigmaPpp.getPhiAtS(s) / pow2(s) / SI::mbarn << "\n";
    }
  }
  {
    auto eAxis = utils::LogAxis(0.1 * SI::GeV, 1e11 * SI::GeV, 10000);
    utils::OutputFile out("test_xsecs_epsprime.txt");
    out << "#\n";
    out << std::scientific;
    for (auto e : eAxis) {
      out << e / SI::GeV << "\t";
      out << sigmaPpp.getAtEpsPrime(e) / SI::mbarn << "\n";
    }
  }
}

void plot_photodisintegration() {
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
}

int main() {
  try {
    utils::startup_information();
    plot_pionproduction();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}