#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H

#include <string>
#include <vector>

#include "simprop/photonFields/PhotonField.h"
#include "simprop/units.h"

namespace simprop {
namespace photonfields {

class Dominguez2011PhotonField final : public PhotonField {
 protected:
  size_t m_zSize;
  size_t m_eSize;
  std::string m_filename;
  std::vector<double> m_redshifts;
  std::vector<double> m_photonEnergies;
  std::vector<double> m_density;
  std::vector<double> m_Igamma;
  //   utils::LookupTable<18, 50> m_field{"data/EBL_Dominguez2011.txt"};
  //   utils::LookupTable<18, 50> m_Igamma{"data/EBL_Igamma_Dominguez2011.txt"};
  //   const double m_IgammaUnits = 1. / utils::pow<2>(SI::eV) / utils::pow<2>(SI::m3);
  const double m_densityUnits = 1. / SI::eV / SI::m3;

 public:
  Dominguez2011PhotonField(size_t zSize, size_t eSize, std::string filename);
  Dominguez2011PhotonField();
  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
  double getMinPhotonEnergy(double z = 0) const override { return m_ePhotonMin; }
  double getMaxPhotonEnergy(double z = 0) const override { return m_ePhotonMax; }

 protected:
  void loadDataFile();
};

}  // namespace photonfields
}  // namespace simprop

#endif