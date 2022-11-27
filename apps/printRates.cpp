#include <gsl/gsl_statistics.h>

#include "simprop.h"

using namespace simprop;

void plot_rates() {
  const auto Gammas = utils::LogAxis<double>(1e8, 1e16, 8 * 32);
  const auto units = 1. / SI::Mpc;
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  const auto pppebl = std::make_shared<interactions::PhotoPionProduction>(sigma, ebl);

  utils::OutputFile out("test_photopion_rates.txt");
  out << std::scientific;
  for (auto Gamma : Gammas) {
    out << Gamma << "\t";
    out << pppcmb->rate(proton, Gamma) / SI::cLight / units << "\t";
    out << pppebl->rate(proton, Gamma) / SI::cLight / units << "\t";
    out << pppebl->computeRateComoving(Gamma, 1.0) / SI::cLight / units << "\t";
    out << pppebl->computeRateComoving(Gamma, 2.0) / SI::cLight / units << "\t";
    out << "\n";
  }
}

void plot_sampled_s() {
  utils::Timer timer("main timer");
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  RandomNumberGenerator rng = utils::RNG<double>(1234);

  const auto pppcmb = interactions::PhotoPionProduction(sigma, cmb);
  utils::OutputFile out("test_photopion_sample_s.txt");
  for (size_t i = 0; i < 300000; ++i) {
    if (i % 1000 == 0) std::cout << i << "\n";
    out << i << "\t";
    out << pppcmb.sampleS(rng(), 4. * SI::GeV2) / SI::GeV2 << "\t";
    out << pppcmb.sampleS(rng(), 40. * SI::GeV2) / SI::GeV2 << "\t";
    out << "\n";
  }
}

void plot_sampled_epsilon() {
  utils::Timer timer("main timer");
  const auto cmb = std::make_shared<photonfields::CMB>();
  LOGI << "CMB min photon energy = " << cmb->getMinPhotonEnergy() / SI::eV;
  LOGI << "CMB max photon energy = " << cmb->getMaxPhotonEnergy() / SI::eV;
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  LOGI << "EBL min photon energy = " << ebl->getMinPhotonEnergy() / SI::eV;
  LOGI << "EBL max photon energy = " << ebl->getMaxPhotonEnergy() / SI::eV;
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  const auto pppebl = std::make_shared<interactions::PhotoPionProduction>(sigma, ebl);
  {
    const auto photonEnergies = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 10000);
    utils::OutputFile out("test_photopion_epsilon_pdf.txt");
    out << "# eps [eV] - pdf\n";
    out << std::scientific;
    for (auto eps : photonEnergies) {
      out << eps / SI::eV << "\t";
      const double E = 1e19 * SI::eV;
      auto s_max = pow2(SI::protonMassC2) + 4. * E * eps;
      out << cmb->density(eps) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      out << ebl->density(eps) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      double z = 1.0;
      out << cmb->density(eps / (1. + z)) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      out << ebl->density(eps / (1. + z), z) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      z = 2.0;
      out << cmb->density(eps / (1. + z)) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      out << ebl->density(eps / (1. + z), z) / pow2(eps) * sigma->getPhiAtS(s_max) << "\t";
      out << "\n";
    }
  }
  {
    RandomNumberGenerator rng = utils::RNG<double>(1234);
    utils::OutputFile out("test_photopion_epsilon_sample.txt");
    for (size_t i = 0; i < 100000; ++i) {
      if (i % 1000 == 0) std::cout << i << "\n";
      out << i << "\t";
      const double E = 1e19 * SI::eV;
      out << pppcmb->sampleEps(rng(), E, 0.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), E, 0.) / SI::eV << "\t";
      out << pppcmb->sampleEps(rng(), E, 1.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), E, 1.) / SI::eV << "\t";
      out << pppcmb->sampleEps(rng(), E, 2.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), E, 2.) / SI::eV << "\t";
      out << "\n";
    }
  }
}

void plot_pion_energies() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  utils::OutputFile out("test_photopion_pionenergy_sample.txt");
  const size_t N = 10000;
  const double z = 0;
  for (size_t i = 0; i < N; ++i) {
    {
      const auto E = 1e18 * SI::eV;
      const auto Gamma = E / SI::protonMassC2;
      auto finalState = pppcmb->finalState({proton, z, Gamma}, z, rng);
      out << finalState[1].getGamma() * SI::pionMassC2 / E << "\t";
    }
    {
      const auto E = 1e20 * SI::eV;
      const auto Gamma = E / SI::protonMassC2;
      auto finalState = pppcmb->finalState({proton, z, Gamma}, z, rng);
      out << finalState[1].getGamma() * SI::pionMassC2 / E << "\n";
    }
  }
}

void plot_inelasticity() {
  RandomNumberGenerator rng = utils::RNG<double>(234);
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, cmb);
  utils::OutputFile out("test_photopion_inelasticity.txt");
  const auto protonEnergies = utils::LogAxis<double>(1e18 * SI::eV, 1e21 * SI::eV, 3 * 8);
  const size_t N = 10000;
  const double z = 0;
  for (auto& E : protonEnergies) {
    double Y[N];
    for (size_t i = 0; i < N; ++i) {
      auto finalState = pppcmb->finalState({proton, z, E / SI::protonMassC2}, z, rng);
      Y[i] = finalState[1].getGamma() * SI::pionMassC2 / E;
    }
    out << E / SI::eV << " ";
    out << gsl_stats_mean(Y, 1, N) << " " << gsl_stats_sd(Y, 1, N) << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    plot_rates();
    plot_sampled_s();
    plot_sampled_epsilon();
    plot_pion_energies();
    plot_inelasticity();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
