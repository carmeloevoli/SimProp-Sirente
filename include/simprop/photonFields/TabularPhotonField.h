#ifndef SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simprop/Units.h"
#include "simprop/utils/interpolators.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace photonfield {

template <size_t zSize, size_t eSize>
class TabularPhotonField {
 public:
  explicit TabularPhotonField(std::string filePath) : m_filePath(filePath) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for photon field does not exist");
    readRedshift();
    readPhotonEnergy();
    readPhotonDensity();
  }

  double get(double ePhoton, double z = 0.) const {
    double logEnergy = std::log10(ePhoton / SI::eV);
    return utils::interpolate2d(z, logEnergy, m_redshifts, m_photonEnergies, m_photonDensity);
  }

 protected:
  void readRedshift() {
    auto v = utils::loadRow(m_filePath, 0, ",");
    if (v.size() != zSize) throw std::runtime_error("error in reading redshifts");
    m_redshifts.reserve(zSize);
    std::copy(v.begin(), v.end(), std::back_inserter(m_redshifts));
  }
  void readPhotonEnergy() {
    auto v = utils::loadRow(m_filePath, 1, ",");
    if (v.size() != eSize) throw std::runtime_error("error in reading photon energies");
    m_photonEnergies.reserve(eSize);
    std::copy(v.begin(), v.end(), std::back_inserter(m_photonEnergies));
  }
  void readPhotonDensity() {
    m_photonEnergies.reserve(eSize * zSize);
    for (size_t i = 0; i < zSize; ++i) {
      auto v = utils::loadRow(m_filePath, i + 2, ",");
      if (v.size() != eSize) throw std::runtime_error("error in reading photon energies");
      std::copy(v.begin(), v.end(), std::back_inserter(m_photonDensity));
    }
  }

 protected:
  std::vector<double> m_redshifts;
  std::vector<double> m_photonEnergies;
  std::vector<double> m_photonDensity;
  std::string m_filePath;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_TABULARPHOTONFIELD_H