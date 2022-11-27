#include "simprop/photonFields/LookupTablePhotonField.h"

#include "simprop/core/common.h"
#include "simprop/utils/io.h"
#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace photonfields {

LookupTablePhotonField::LookupTablePhotonField(size_t zSize, size_t eSize, std::string filename) {
  m_zSize = zSize;
  m_eSize = eSize;
  m_filename = "data/" + filename;
  auto fileSize = utils::countFileLines(m_filename);
  if (utils::fileExists(m_filename) && fileSize == (zSize * eSize)) {
    m_redshifts.reserve(zSize);
    m_logPhotonEnergies.reserve(eSize);
    m_logDensity.reserve(zSize * eSize);
    m_logIgamma.reserve(zSize * eSize);
    loadDataFile();
  } else {
    throw std::runtime_error("error reading from file : " + filename);
  }
  LOGD << "calling " << __func__ << " constructor";
}

void LookupTablePhotonField::loadDataFile() {
  auto v = utils::loadFileByRow(m_filename, ",");
  size_t counter = 0;
  for (size_t i = 0; i < m_zSize; ++i) {
    for (size_t j = 0; j < m_eSize; ++j) {
      auto line = v.at(counter);
      if (line.size() != 4) throw std::runtime_error("error in reading table values");
      if (j == 0) m_redshifts.emplace_back(line[0]);
      if (i == 0) m_logPhotonEnergies.emplace_back(line[1]);
      m_logDensity.emplace_back(line[2]);
      m_logIgamma.emplace_back(line[3]);
      counter++;
    }
  }
  assert(m_redshifts.size() == m_zSize && m_logPhotonEnergies.size() == m_eSize);
}

double LookupTablePhotonField::density(double ePhoton, double z) const {
  constexpr auto units = 1. / SI::eV / SI::m3;
  double value = 0;
  auto loge = std::log10(ePhoton / SI::eV);
  if (utils::isInside(z, m_redshifts) && utils::isInside(loge, m_logPhotonEnergies)) {
    auto logn = utils::interpolate2d(z, loge, m_redshifts, m_logPhotonEnergies,
                                     m_logDensity);  // TODO speed up this
    value = std::pow(10., logn) * units;
  }
  return std::max(value, 0.);
}

double LookupTablePhotonField::I_gamma(double ePhoton, double z) const {
  constexpr auto units = 1. / pow2(SI::eV) / SI::m3;
  double value = 0;
  auto loge = std::log10(ePhoton / SI::eV);
  if (utils::isInside(z, m_redshifts) && utils::isInside(loge, m_logPhotonEnergies)) {
    auto logn = utils::interpolate2d(z, loge, m_redshifts, m_logPhotonEnergies, m_logIgamma);
    value = std::pow(10., logn) * units;
  }
  return value;
}

}  // namespace photonfields
}  // namespace simprop