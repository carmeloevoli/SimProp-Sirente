// Copyright 2023 SimProp-dev [MIT License]
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
  using std::log10;
  using std::max;
  auto v = utils::loadFileByRow(m_filename, ",");
  size_t counter = 0;
  for (size_t i = 0; i < m_zSize; ++i) {
    for (size_t j = 0; j < m_eSize; ++j) {
      auto line = v.at(counter);
      if (line.size() != 4) throw std::runtime_error("error in reading table values");
      auto z = line[0];
      auto eps = line[1] * SI::eV;
      auto n = line[2] * (1. / SI::eV / SI::m3);
      auto I_gamma = line[3] * (1. / pow2(SI::eV) / SI::m3);
      if (j == 0) m_redshifts.emplace_back(z);
      if (i == 0) m_logPhotonEnergies.emplace_back(log10(eps));
      m_logDensity.emplace_back(log10(max(n, 1e-30)));
      m_logIgamma.emplace_back(log10(max(I_gamma, 1e-30)));
      counter++;
    }
  }
  assert(m_redshifts.size() == m_zSize && m_logPhotonEnergies.size() == m_eSize);
}

double LookupTablePhotonField::density(double epsRestFrame, double z) const {
  double value = 0;
  auto loge = std::log10(epsRestFrame);
  if (utils::isInside(z, m_redshifts) && utils::isInside(loge, m_logPhotonEnergies)) {
    auto logn = utils::interpolate2d(z, loge, m_redshifts, m_logPhotonEnergies,
                                     m_logDensity);  // TODO(CE) speed up this
    value = std::pow(10., logn);
  }
  return std::max(value, 0.);
}

double LookupTablePhotonField::I_gamma(double epsRestFrame, double z) const {
  double value = 0;
  auto loge = std::log10(epsRestFrame);
  if (utils::isInside(z, m_redshifts) && utils::isInside(loge, m_logPhotonEnergies)) {
    auto logn = utils::interpolate2d(z, loge, m_redshifts, m_logPhotonEnergies, m_logIgamma);
    value = std::pow(10., logn);
  }
  return value;
}

}  // namespace photonfields
}  // namespace simprop
