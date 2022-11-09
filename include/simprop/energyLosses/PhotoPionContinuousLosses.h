#ifndef SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H
#define SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H

#include "simprop/crossSections/PhotoPionProductionXsec.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace losses {

double inelasticity(double s);
double inelasticityPoorApproximation(double s, double cosTheta_pi);

class PhotoPionContinuousLosses final : public ContinuousLosses {
 protected:
  photonfields::PhotonFields m_photonFields;
  std::shared_ptr<xsecs::PhotoPionProductionXsec> m_sigma;

 public:
  PhotoPionContinuousLosses(const std::shared_ptr<photonfields::PhotonField>& photonField);
  PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields);
  virtual ~PhotoPionContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double computeBetaComoving(double Gamma) const;
};

}  // namespace losses

}  // namespace simprop

#endif