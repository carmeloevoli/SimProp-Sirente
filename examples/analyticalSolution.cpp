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

  double energyJacobian(double E, double zMax, double relError = 1e-3) const {
    auto dydz = [&](double z, double y) {
      auto E_g = sourceEnergy(E, z, 1e-5);
      auto E_prime = E_g * (1. + z);
      auto beta = m_pp->betaComoving(E_prime / SI::protonMassC2);
      auto dbeta = utils::deriv<double>(
          [&](double x) { return m_pp->betaComoving(x / SI::protonMassC2); }, E_prime, 1e-4);
      return y / (1. + z) + y * m_cosmology->dtdz(z) * pow3(1. + z) * (beta + E_prime * dbeta);
    };
    return utils::odeiv<double>(dydz, 1., 0., zMax, relError);
  }
};

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer for analytical solution");
    Berezinsky b;
    auto z = utils::LogAxis(1e-3, 1., 100);
    {
      utils::OutputFile out("proton_characteristics.txt");
      for (const auto &z_i : z) {
        std::cout << z_i << "\n";
        out << std::scientific << z_i << "\t";
        out << b.sourceEnergy(1e17 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
        out << b.sourceEnergy(1e18 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
        out << b.sourceEnergy(1e19 * SI::eV, z_i, 1e-4) / SI::eV << "\t";
        out << "\n";
      }
    }
    {
      utils::OutputFile out("proton_jacobian.txt");
      for (const auto &z_i : z) {
        std::cout << z_i << "\n";
        out << std::scientific << z_i << "\t";
        out << b.energyJacobian(1e17 * SI::eV, z_i, 1e-3) << "\t";
        out << b.energyJacobian(1e18 * SI::eV, z_i, 1e-3) << "\t";
        out << b.energyJacobian(1e19 * SI::eV, z_i, 1e-3) << "\t";
        out << "\n";
      }
    }
  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
