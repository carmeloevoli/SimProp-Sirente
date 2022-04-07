#include <gsl/gsl_statistics.h>

#include "simprop.h"

using namespace simprop;

void plot_rates() {
  const auto Gammas = utils::LogAxis(1e8, 1e16, 8 * 32);
  const auto units = 1. / SI::Mpc;
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  const auto pppebl = std::make_shared<interactions::PhotoPionProduction>(sigma, ebl);

  utils::OutputFile out("test_rates.txt");
  out << std::scientific;
  for (auto Gamma : Gammas) {
    out << Gamma << "\t";
    out << pppcmb->rate(proton, Gamma) / SI::cLight / units << "\t";
    out << pppebl->rate(proton, Gamma) / SI::cLight / units << "\t";
    out << pppebl->rateCrpropaComparison(proton, Gamma, 1.0) / SI::cLight / units << "\t";
    out << pppebl->rateCrpropaComparison(proton, Gamma, 2.0) / SI::cLight / units << "\t";
    out << "\n";
  }
}

void plot_sampled_s() {
  // utils::Timer timer("main timer");
  // const auto cmb = std::make_shared<photonfields::CMB>();
  // const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  // RandomNumberGenerator rng = utils::RNG<double>(1234);

  // const auto pppcmb = interactions::PhotoPionProduction(sigma, cmb, rng);
  // utils::OutputFile out("test_sample_s.txt");
  // for (size_t i = 0; i < 300000; ++i) {
  //   if (i % 1000 == 0) std::cout << i << "\n";
  //   out << i << "\t";
  //   out << pppcmb.sample_s(rng(), 4. * SI::GeV2) / SI::GeV2 << "\t";
  //   out << pppcmb.sample_s(rng(), 40. * SI::GeV2) / SI::GeV2 << "\t";
  //   out << "\n";
  // }
}

void plot_epsilon_pdf() {
  // utils::Timer timer("main timer");
  // const auto cmb = std::make_shared<photonfields::CMB>();
  // LOGI << "min photon energy = " << cmb->getMinPhotonEnergy() / SI::eV;
  // LOGI << "max photon energy = " << cmb->getMaxPhotonEnergy() / SI::eV;
  // const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  // const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  // {
  //   const auto photonEnergies = utils::LogAxis(1e-5 * SI::eV, 1e-2 * SI::eV, 2000);
  //   utils::OutputFile out("test_epsilon_pdf.txt");
  //   out << "# eps [eV] - pdf\n";
  //   out << std::scientific;
  //   for (auto eps : photonEnergies) {
  //     out << eps / SI::eV << "\t";
  //     auto s_max = pow2(SI::protonMassC2) + 4. * (1e22 * SI::eV) * eps;
  //     out << cmb->density(eps) / pow2(eps) * pppcmb->phi(s_max) << "\t";
  //     out << "\n";
  //   }
  // }
  // {
  //   RandomNumberGenerator rng = utils::RNG<double>(1234);
  //   utils::OutputFile out("test_epsilon_sample.txt");
  //   for (size_t i = 0; i < 3000; ++i) {
  //     if (i % 100 == 0) std::cout << i << "\n";
  //     out << i << "\t";
  //     out << pppcmb->sample_eps(rng(), 1e22 * SI::eV, 0.) / SI::eV << "\t";
  //     // out << pppcmb.sample_s(rng(), 40. * SI::GeV2) / SI::GeV2 << "\t";
  //     out << "\n";
  //   }
  // }
}

void plot_pion_energies() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  utils::OutputFile out("test_pionenergy_sample.txt");
  const size_t N = 1000;
  const double E18 = 1e18 * SI::eV;
  const double E20 = 1e20 * SI::eV;
  for (size_t i = 0; i < N; ++i) {
    out << pppcmb->samplePionEnergy(E18, 0., rng) / E18 << " ";
    out << pppcmb->samplePionEnergy(E20, 0., rng) / E20 << "\n";
  }
}

void plot_inelasticity() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  utils::OutputFile out("test_inelasticity.txt");
  const auto protonEnergies = utils::LogAxis(1e18 * SI::eV, 1e21 * SI::eV, 4 * 8);
  const size_t N = 100;
  for (auto& E : protonEnergies) {
    double Y[N];
    for (size_t i = 0; i < N; ++i) {
      Y[i] = pppcmb->samplePionEnergy(E, 0., rng) / E;
    }
    out << E / SI::eV << " ";
    out << gsl_stats_median(Y, 1, N) << " " << gsl_stats_sd(Y, 1, N) << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    // plot_rates();
    // plot_sampled_s();
    // plot_epsilon_pdf();
    plot_pion_energies();
    plot_inelasticity();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
