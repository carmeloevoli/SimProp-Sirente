#include "simprop.h"

using namespace simprop;

void plot_pionproduction() {
  auto sigmaPhotoPion = xsecs::PhotoPionXsec();
  auto epsPrimeAxis = utils::LogAxis<double>(0.1 * SI::GeV, 1e11 * SI::GeV, 10000);
  {
    utils::OutputFile out("test_photopion_xsecs.txt");
    out << "# eps [GeV] - s [GeV^2] - sigma proton [mbarn] - sigma neutron [mbarn]\n";
    out << std::scientific;
    for (auto epsPrime : epsPrimeAxis) {
      auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
      out << epsPrime / SI::GeV << "\t";
      out << s / SI::GeV2 << "\t";
      out << sigmaPhotoPion.getAtS(proton, s) / SI::mbarn << "\t";
      out << sigmaPhotoPion.getAtS(neutron, s) / SI::mbarn << "\t";
      out << "\n";
    }
  }
  {
    utils::OutputFile out("test_photopion_phi.txt");
    out << "#\n";
    out << std::scientific;
    for (auto epsPrime : epsPrimeAxis) {
      auto s = pow2(SI::protonMassC2) + 2. * SI::protonMassC2 * epsPrime;
      out << epsPrime / SI::GeV << "\t";
      out << s / SI::GeV2 << "\t";
      out << sigmaPhotoPion.getPhiAtS(proton, s) / pow2(s) / SI::mbarn << "\t";
      out << sigmaPhotoPion.getPhiAtS(neutron, s) / pow2(s) / SI::mbarn << "\t";
      out << "\n";
    }
  }
}

void plot_photodisintegration() {
  // {
  //   auto sigmaPd = interactions::PsbPhotoDisintegration();
  //   auto eAxis = utils::LinAxis<double>(0.1 * SI::MeV, 1e2 * SI::MeV, 1000);
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