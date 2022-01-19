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
  std::vector<double> m_logPhotonEnergies;
  std::vector<double> m_logDensity;
  std::vector<double> m_logIgamma;

 public:
  Dominguez2011PhotonField(size_t zSize, size_t eSize, std::string filename);
  Dominguez2011PhotonField();
  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
  double getMinPhotonEnergy(double z = 0) const override {
    return std::pow(10., m_logPhotonEnergies.front()) * SI::eV;
  }
  double getMaxPhotonEnergy(double z = 0) const override {
    return std::pow(10., m_logPhotonEnergies.back()) * SI::eV;
  }

 protected:
  void loadDataFile();
};

}  // namespace photonfields
}  // namespace simprop

#endif