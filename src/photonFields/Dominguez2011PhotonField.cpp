#include "simprop/photonFields/Dominguez2011PhotonField.h"

#include "simprop/common.h"
#include "simprop/utils/io.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace photonfield {

Dominguez2011PhotonField::Dominguez2011PhotonField(size_t zSize, size_t eSize, std::string filename)
    : m_zSize(zSize), m_eSize(eSize), m_filename("data/" + filename) {
  auto fileSize = utils::countFileLines(m_filename);
  if (utils::fileExists(m_filename) && fileSize == (zSize * eSize)) {
    m_redshifts.reserve(zSize);
    m_photonEnergies.reserve(eSize);
    m_density.reserve(zSize * eSize);
    m_Igamma.reserve(zSize * eSize);
    loadDataFile();
  } else {
    throw std::runtime_error("error reading from file : " + filename);
  }
}

Dominguez2011PhotonField::Dominguez2011PhotonField()
    : Dominguez2011PhotonField(18, 50, "EBL_Dominguez2011.txt") {}

void Dominguez2011PhotonField::loadDataFile() {
  auto v = utils::loadFileByRow(m_filename, ",");
  size_t counter = 0;
  for (size_t i = 0; i < m_zSize; ++i) {
    for (size_t j = 0; j < m_eSize; ++j) {
      auto line = v.at(counter);
      if (line.size() != 3) throw std::runtime_error("error in reading table values");
      if (j == 0) m_redshifts.emplace_back(line[0]);
      if (i == 0) m_photonEnergies.emplace_back(line[1]);
      m_density.emplace_back(line[2]);
      m_Igamma.emplace_back(0.);  // TODO this
      counter++;
    }
  }
  assert(m_redshifts.size() == m_zSize && m_photonEnergies.size() == m_eSize);
}

double Dominguez2011PhotonField::density(double ePhoton, double z) const {
  double value = 0;
  auto loge = std::log10(ePhoton / SI::eV);
  if (utils::isInside(z, m_redshifts) && utils::isInside(loge, m_photonEnergies)) {
    auto logDensity = utils::interpolate2d(z, loge, m_redshifts, m_photonEnergies, m_density);
    value = std::pow(10., logDensity) * m_densityUnits;
  }
  return value;
};

double Dominguez2011PhotonField::I_gamma(double ePhoton, double z) const {
  // const auto loge = std::log10(ePhoton / SI::eV);
  // if (m_field.isWithinXRange(z) && m_field.isWithinYRange(loge)) {
  //   const auto value = m_Igamma.get(z, loge);
  //   const auto I = std::pow(10., value) * m_IgammaUnits;
  //   return I;
  // } else {
  //   return 0;
  // }
  return 0;
}

}  // namespace photonfield
}  // namespace simprop