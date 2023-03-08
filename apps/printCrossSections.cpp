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
  auto sigma = xsecs::PhotoDisintegrationTalysXsec();
  auto epsAxis = utils::LinAxis<double>(0.1 * SI::MeV, 1e2 * SI::MeV, 1000);
  utils::OutputFile out("test_pd_v2r4_sigma.txt");
  out << std::scientific;
  for (auto eps : epsAxis) {
    out << eps / SI::MeV << "\t";
    out << sigma.getSingleNucleon(Fe56, eps) / SI::mbarn << "\t";
    out << sigma.getSingleNucleon(O16, eps) / SI::mbarn << "\t";
    out << sigma.getSingleNucleon(C12, eps) / SI::mbarn << "\t";
    out << sigma.getAlpha(Fe56, eps) / SI::mbarn << "\t";
    out << sigma.getAlpha(O16, eps) / SI::mbarn << "\t";
    out << sigma.getAlpha(C12, eps) / SI::mbarn << "\t";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    // plot_pionproduction();
    plot_photodisintegration();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}