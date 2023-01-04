#ifndef SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H
#define SIMPROP_LOSSES_PHOTOPION_CONTINUOUS_H

#include "simprop/crossSections/PhotoPionXsecs.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace losses {

double inelasticity(double s);
double inelasticityPoorApproximation(double s, double cosTheta_pi);

class PhotoPionContinuousLosses final : public ContinuousLosses {
 protected:
  photonfields::PhotonFields m_photonFields;
  xsecs::PhotoPionProtonXsec m_xs_proton;
  xsecs::PhotoPionNeutronXsec m_xs_neutron;

 public:
  PhotoPionContinuousLosses(const std::shared_ptr<photonfields::PhotonField>& photonField);
  PhotoPionContinuousLosses(const photonfields::PhotonFields& photonFields);
  virtual ~PhotoPionContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;

 protected:
  double computeBetaComoving(double Gamma, double z, const xsecs::CrossSection& xs) const;
};

}  // namespace losses

}  // namespace simprop

#endif