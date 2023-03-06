// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_ENERGYLOSSES_PHOTOPIONCONTINUOUSLOSSES_H_
#define SIMPROP_ENERGYLOSSES_PHOTOPIONCONTINUOUSLOSSES_H_

#include <memory>

#include "simprop/crossSections/PhotoPionXsecs.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/PhotonField.h"

namespace simprop {
namespace losses {

enum class MODE { ISOTROPIC, BACKWARD, SOPHIA };

double inelasticity(double s, MODE m);

class PhotoPionContinuousLosses final : public ContinuousLosses {
 protected:
  photonfields::PhotonFields m_photonFields;
  xsecs::PhotoPionXsec m_xs;

 public:
  PhotoPionContinuousLosses(const std::shared_ptr<photonfields::PhotonField>& photonField,
                            MODE meanAngle = MODE::ISOTROPIC);
  PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields,
                            MODE meanAngle = MODE::ISOTROPIC);
  virtual ~PhotoPionContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  MODE m_meanAngle;
};

}  // namespace losses
}  // namespace simprop

#endif  // SIMPROP_ENERGYLOSSES_PHOTOPIONCONTINUOUSLOSSES_H_
