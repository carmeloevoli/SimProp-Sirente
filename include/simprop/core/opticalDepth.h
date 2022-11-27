#ifndef SIMPROP_CORE_OPTICALDEPTH_H
#define SIMPROP_CORE_OPTICALDEPTH_H

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace core {

class OpticalDepth {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_ebl;
  std::shared_ptr<photonfields::CMB> m_cmb;

 public:
  OpticalDepth(const std::shared_ptr<cosmo::Cosmology>& cosmology,
               const std::shared_ptr<photonfields::PhotonField>& phField)
      : m_cosmology(cosmology), m_ebl(phField) {
    m_cmb = std::make_shared<photonfields::CMB>();
  }

  double get(double eGamma, double zSource) const;

  virtual ~OpticalDepth() = default;

 protected:
  double integrateOverAngle(double eGamma, double z) const;
  double integrateOverPhField(double eGamma, double z, double mu) const;
};

}  // namespace core
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
