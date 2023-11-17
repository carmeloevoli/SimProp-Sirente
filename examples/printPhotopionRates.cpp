#include <gsl/gsl_statistics.h>

#include "simprop.h"

using namespace simprop;

void plot_rates() {
  const auto Gammas = utils::LogAxis<double>(1e8, 1e16, 8 * 32);
  const auto units = 1. / SI::Mpc;
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto ebl = std::make_shared<photonfields::Saldana2021PhotonField>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(cmb);
  const auto pppebl = std::make_shared<interactions::PhotoPionProduction>(ebl);
  {
    utils::OutputFile out("SimProp_photopion_rates.txt");
    out << std::scientific;
    for (auto Gamma : Gammas) {
      out << Gamma << "\t";
      out << pppcmb->rate(proton, Gamma) / SI::cLight / units << "\t";
      out << pppebl->rate(proton, Gamma) / SI::cLight / units << "\t";
      out << pppcmb->rate(proton, Gamma, 1.) / SI::cLight / units << "\t";
      out << pppebl->rate(proton, Gamma, 1.) / SI::cLight / units << "\t";
      out << pppcmb->rate(proton, Gamma, 3.) / SI::cLight / units << "\t";
      out << pppebl->rate(proton, Gamma, 3.) / SI::cLight / units << "\t";
      out << pppcmb->rate(proton, Gamma, 5.) / SI::cLight / units << "\t";
      out << pppebl->rate(proton, Gamma, 5.) / SI::cLight / units << "\t";
      out << "\n";
    }
  }
}

void plot_sampled_s() {
  utils::Timer timer("sampling s");
  const auto cmb = std::make_shared<photonfields::CMB>();
  const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(cmb);

  RandomNumberGenerator rng = utils::RNG<double>(1234);
  utils::OutputFile out("SimProp_photopion_sample_s.txt");
  for (size_t i = 0; i < 300000; ++i) {
    if (i % 1000 == 0) std::cout << i << "\n";
    out << i << "\t";
    out << pppcmb->sampleS(rng(), proton, 4. * SI::GeV2) / SI::GeV2 << "\t";
    out << pppcmb->sampleS(rng(), proton, 40. * SI::GeV2) / SI::GeV2 << "\t";
    out << "\n";
  }
}

void plot_sampled_epsilon() {
  utils::Timer timer("sampling epsilon");
  const auto cmb = std::make_shared<photonfields::CMB>();
  LOGI << "CMB min photon energy = " << cmb->getMinPhotonEnergy() / SI::eV;
  LOGI << "CMB max photon energy = " << cmb->getMaxPhotonEnergy() / SI::eV;
  const auto ebl = std::make_shared<photonfields::Dominguez2011PhotonField>();
  LOGI << "EBL min photon energy = " << ebl->getMinPhotonEnergy() / SI::eV;
  LOGI << "EBL max photon energy = " << ebl->getMaxPhotonEnergy() / SI::eV;
  const auto xs = std::make_shared<xsecs::PhotoPionXsec>();
  {
    const auto photonEnergies = utils::LogAxis<double>(1e-3 * SI::eV, 1e1 * SI::eV, 1000);
    std::vector<double> redshifts = {0., 1., 2., 3.};
    utils::OutputFile out("SimProp_photopion_epsilon_pdf.txt");
    out << "# eps [eV] - pdf\n";
    out << std::scientific;
    for (auto eps : photonEnergies) {
      out << eps / SI::eV << "\t";
      const double E = 1e19 * SI::eV;
      auto s_max = pow2(SI::protonMassC2) + 4. * E * eps;
      for (auto z : redshifts) {
        out << cmb->density(eps, z) / pow2(eps) * xs->getPhiAtS(proton, s_max) << "\t";
        out << ebl->density(eps, z) / pow2(eps) * xs->getPhiAtS(proton, s_max) << "\t";
      }
      out << "\n";
    }
  }
  {
    RandomNumberGenerator rng = utils::RNG<double>(138);
    const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(cmb);
    const auto pppebl = std::make_shared<interactions::PhotoPionProduction>(ebl);

    utils::OutputFile out("SimProp_photopion_epsilon_sample.txt");
    for (size_t i = 0; i < 10000; ++i) {
      if (i % 1000 == 0) {
        LOGD << i;
      }
      out << i << "\t";
      const double E = 1e19 * SI::eV;
      out << pppcmb->sampleEps(rng(), proton, E, 0.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), proton, E, 0.) / SI::eV << "\t";
      out << pppcmb->sampleEps(rng(), proton, E, 1.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), proton, E, 1.) / SI::eV << "\t";
      out << pppcmb->sampleEps(rng(), proton, E, 2.) / SI::eV << "\t";
      out << pppebl->sampleEps(rng(), proton, E, 2.) / SI::eV << "\t";
      out << "\n";
    }
  }
}

// double sum_pion_energy(const ParticleStack& stack) {
//   double pion_energy = 0.;
//   for (const auto& particle : stack) {
//     if (pidIsPion(particle.getPid())) pion_energy += particle.getGamma() * SI::pionMassC2;
//   }

auto getPionsEnergy = [](double E, const Particle& p) {
  auto value = (pidIsPion(p.getPid())) ? p.getGamma() * SI::pionMassC2 : 0.;
  return E + value;
};

void plot_pion_energies() {
  const auto cmb = std::make_shared<photonfields::CMB>();
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  const size_t N = 100000;
  const double z = 0;
  std::vector<double> protonEnergies = {1e18 * SI::eV, 1e19 * SI::eV, 1e20 * SI::eV, 1e21 * SI::eV};
  {
    const auto pppcmb = std::make_shared<interactions::PhotoPionProduction>(cmb);
    utils::OutputFile out("SimProp_photopion_pionenergy_sample.txt");
    for (size_t i = 0; i < N; ++i) {
      for (auto E : protonEnergies) {
        const auto Gamma = E / SI::protonMassC2;
        auto finalState = pppcmb->finalState({proton, z, Gamma}, z, rng);
        double piEnergy = std::accumulate(finalState.begin(), finalState.end(), 0., getPionsEnergy);
        out << piEnergy / E << "\t";
      }
      out << "\n";
    }
  }
  {
    const auto pppcmb = std::make_shared<interactions::PhotoPionProductionSophia>(cmb);
    utils::OutputFile out("SimProp_photopionsophia_pionenergy_sample.txt");
    for (size_t i = 0; i < N; ++i) {
      for (auto E : protonEnergies) {
        const auto Gamma = E / SI::protonMassC2;
        auto finalState = pppcmb->finalState({proton, z, Gamma}, z, rng);
        auto piEnergy = std::accumulate(finalState.begin(), finalState.end(), 0., getPionsEnergy);
        out << piEnergy / E << "\t";
      }
      out << "\n";
    }
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    plot_rates();
    plot_sampled_s();
    plot_sampled_epsilon();
    plot_pion_energies();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
