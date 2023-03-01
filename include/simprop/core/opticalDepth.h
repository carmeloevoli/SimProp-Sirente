// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_CORE_OPTICALDEPTH_H_
#define SIMPROP_CORE_OPTICALDEPTH_H_

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace core {

class OpticalDepth {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  photonfields::PhotonFields m_photonFields;

 public:
  OpticalDepth(const std::shared_ptr<cosmo::Cosmology>& cosmology,
               const std::shared_ptr<photonfields::PhotonField>& photonField);
  OpticalDepth(const std::shared_ptr<cosmo::Cosmology>& cosmology,
               const photonfields::PhotonFields& photonFields);

  double get(double eGamma, double zSource) const;

  virtual ~OpticalDepth() = default;

 protected:
  double integrateOverAngle(double eGamma, double z) const;
  double integrateOverPhField(double eGamma, double z, double mu) const;
};

}  // namespace core
}  // namespace simprop

#endif  // SIMPROP_CORE_OPTICALDEPTH_H_
