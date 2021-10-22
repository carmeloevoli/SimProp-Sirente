#ifndef SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H

#include <array>
#include <stdexcept>

#include "simprop/utils/misc.h"

namespace simprop {
namespace photonfield {

template <size_t zSize, size_t eSize>
class TabularPhotonField {
 public:
  explicit TabularPhotonField(std::string filePath) : m_filePath(filePath) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for photon field does not exist");
  }

  double get(double ePhoton, double z = 0.) const {
    if (ePhoton < m_photonEnergies.front() || ePhoton > m_photonEnergies.back())
      return 0;
    else if (z < m_redshifts.front() || z > m_redshifts.back())
      return 0;
    else
      return 1;
  }

 protected:
  // void readPhotonEnergy(std::string filePath);
  // void readPhotonDensity(std::string filePath);
  // void readRedshift(std::string filePath);
  // void checkInputData() const;

  std::array<double, zSize> m_redshifts;
  std::array<double, eSize> m_photonEnergies;
  std::array<double, eSize * zSize> m_photonDensity;
  std::string m_filePath;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H