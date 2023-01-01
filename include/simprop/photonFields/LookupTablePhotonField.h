#ifndef SIMPROP_PHOTONFIELDS_LOOKUPTABLE_H
#define SIMPROP_PHOTONFIELDS_LOOKUPTABLE_H

#include <string>
#include <vector>

#include "simprop/core/units.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace photonfields {

class LookupTablePhotonField : public PhotonField {
 protected:
  size_t m_zSize, m_eSize;
  std::string m_filename;
  std::vector<double> m_redshifts;
  std::vector<double> m_logPhotonEnergies;
  std::vector<double> m_logDensity;
  std::vector<double> m_logIgamma;

 public:
  LookupTablePhotonField(size_t zSize, size_t eSize, std::string filename);

  double density(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;

  double getMinPhotonEnergy() const override { return std::pow(10., m_logPhotonEnergies.front()); }
  double getMaxPhotonEnergy() const override { return std::pow(10., m_logPhotonEnergies.back()); }

 protected:
  void loadDataFile();
};

}  // namespace photonfields
}  // namespace simprop

#endif