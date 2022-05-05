#include "simprop.h"

using namespace simprop;

class Berezinsky {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_cmb;
  std::shared_ptr<losses::PairProductionLosses> m_pp;
  std::shared_ptr<losses::AdiabaticContinuousLosses> m_ad;
  // std::vector<std::shared_ptr<losses::ContinuousLosses> > m_continuousLosses;
  //  std::shared_ptr<interactions::PhotoPionProduction> m_pppcmb;

  double m_slope = 2.7;
  double m_sourceEmissivity = 1.;
  double m_sourceEvolution = 3.;

 public:
  Berezinsky() {
    m_cosmology = std::make_shared<cosmo::Planck2018>();
    m_cmb = std::make_shared<photonfields::CMB>();
    m_pp = std::make_shared<losses::PairProductionLosses>(photonfields::PhotonFields{m_cmb});
    m_ad = std::make_shared<losses::AdiabaticContinuousLosses>(m_cosmology);
  };

  virtual ~Berezinsky() = default;

  double sourceEnergy(double E, double zMax, double relError = 1e-3) const {
    auto dydx = [&](double z, double E_g) {
      return E_g / (1. + z) + E_g * m_cosmology->dtdz(z) * pow3(1. + z) *
                                  m_pp->betaComoving(E_g / SI::protonMassC2 * (1. + z));
    };
    return utils::odeiv<double>(dydx, E, 0., zMax, relError);
  }

  double jacobian(double E, double zMax, double relError = 1e-3) const {
    auto dydz = [&](double z, double y) {
      auto E_g = sourceEnergy(E, z, 1e-4);
      auto E_prime = E_g * (1. + z);
      auto beta = m_pp->betaComoving(E_prime / SI::protonMassC2);
      auto dbeta = utils::deriv5pt<double>(
          [&](double x) { return m_pp->betaComoving(x / SI::protonMassC2); }, E_prime,
          1e-2 * E_prime);
      return y / (1. + z) + y * m_cosmology->dtdz(z) * pow3(1. + z) * (beta + E_prime * dbeta);
    };
    return utils::odeiv<double>(dydz, 1., 0., zMax, relError);
  }

  double computeFlux(double E, double zMax) const {
    const auto K = m_slope - 2;
    const auto L_0 = m_sourceEmissivity;
    const auto m = m_sourceEvolution;
    const auto E_0 = 1e17 * SI::eV;
    const auto factor = SI::cLight / 4. / M_PI * K * L_0;
    auto integrand = [&](double z) {
      auto E_g = sourceEnergy(E, z, 1e-4);
      auto y = jacobian(E, z, 1e-4);
      return m_cosmology->dtdz(z) * std::pow(1. + z, m) * std::pow(E_g / E_0, -m_slope) * y;
    };
    auto I = utils::QAGIntegration<double>(integrand, 0., zMax, 1000, 1e-3);
    return factor * I;
  }
};

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for analytical solution");
    Berezinsky b;
    auto z = utils::LogAxis(1e-4, 1.0, 100);
    // {
    //   utils::OutputFile out("proton_characteristics.txt");
    //   for (const auto &z_i : z) {
    //     std::cout << z_i << "\n";
    //     out << std::scientific << z_i << "\t";
    //     out << b.sourceEnergy(1e17 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
    //     out << b.sourceEnergy(1e18 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
    //     out << b.sourceEnergy(1e19 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
    //     out << "\n";
    //   }
    // }
    // {
    //   utils::OutputFile out("proton_jacobian.txt");
    //   for (const auto &z_i : z) {
    //     std::cout << z_i << "\n";
    //     out << std::scientific << z_i << "\t";
    //     out << b.jacobian(1e17 * SI::eV, z_i, 1e-4) << "\t";
    //     out << b.jacobian(1e18 * SI::eV, z_i, 1e-4) << "\t";
    //     out << b.jacobian(1e19 * SI::eV, z_i, 1e-4) << "\t";
    //     out << "\n";
    //   }
    // }
    {
      utils::OutputFile out("proton_spectrum.txt");
      auto E = utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 10);
      for (const auto &E_i : E) {
        std::cout << E_i / SI::eV << "\n";
        out << std::scientific << E_i / SI::eV << "\t";
        out << b.computeFlux(E_i, 1.) << "\t";
        out << "\n";
      }
    }
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
