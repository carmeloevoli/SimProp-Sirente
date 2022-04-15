#ifndef SIMPROP_INTERACTIONS_PPPEBL_H
#define SIMPROP_INTERACTIONS_PPPEBL_H

#include <memory>

#include "simprop/interactions/Interaction.h"
#include "simprop/units.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace interactions {

struct RndUnifNumber {
 public:
  RndUnifNumber(double r) : m_r(r) {
    if (r < 0. || r > 1.) throw std::invalid_argument("invalid random number");
  }
  double get() const { return m_r; };

 protected:
  double m_r;
};

class PhotoPionProduction final : public Interaction {
 protected:
  const double m_sThreshold = pow2(SI::protonMassC2 + SI::pionMassC2);

 public:
  PhotoPionProduction(const std::shared_ptr<xsecs::CrossSection>& sigma,
                      const std::shared_ptr<photonfields::PhotonField>& phField)
      : Interaction(sigma, phField) {}
  virtual ~PhotoPionProduction() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;
  double computeRateComoving(double Gamma, double z) const;  // TODO to make protected
  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;

  double sampleS(RndUnifNumber r, double sMax) const;
  double sampleEps(RndUnifNumber r, double nucleonEnergy, double z) const;
  double samplePionInelasticity(RndUnifNumber r, double s) const;
  double epsPdfIntegral(double photonEnergy, double nucleonEnergy, double z) const;
  PID samplePionCharge(RndUnifNumber r, bool isNeutron) const;

  // protected:
  // double phi(double s) const;
};

}  // namespace interactions
}  // namespace simprop

#endif