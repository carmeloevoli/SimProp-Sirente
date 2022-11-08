#ifndef SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H
#define SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H

#include "simprop/crossSections/PhotoPionProductionXsec.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace losses {

double inelasticity(double s);
double inelasticityPoorApproximation(double s);

class PhotoPionContinuousLosses final : public ContinuousLosses {
 public:
  PhotoPionContinuousLosses();
  virtual ~PhotoPionContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  std::shared_ptr<photonfields::CMB> m_cmb;
  std::shared_ptr<photonfields::PhotonField> m_ebl;
  std::shared_ptr<xsecs::PhotoPionProductionXsec> m_sigma;

 protected:
  double computeBetaComoving(double Gamma, double z) const;
};

}  // namespace losses

}  // namespace simprop

#endif