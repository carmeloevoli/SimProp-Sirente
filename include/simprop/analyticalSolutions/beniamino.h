#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/core/cosmology.h"
#include "simprop/crossSections/KelnerAharonian2008.h"
#include "simprop/energyLosses/PairProductionLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace solutions {

struct BeniaminoParams {
  double slope;
  double sourceEvolution;
  double sourceCutoff;
};

class Beniamino {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  utils::LookupArray<5000> m_losses;

  const double m_sourceEmissivity = 0.4e46 * SI::erg / SI::Mpc3 / SI::year;
  const double m_maxEnergy = 1e23 * SI::eV;
  const double m_minEnergy = 1e17 * SI::eV;

  double m_slope = 2.6;
  double m_sourceEvolution = 0.;
  double m_sourceCutoff = -1.;

 public:
  Beniamino(bool doPhotoPion);
  Beniamino(bool doPhotoPion, BeniaminoParams params);

  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zNow, double zMax, double relError = 1e-3) const;
  double dilationFactor(double E, double zNow, double zMax, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double zMax, double relError = 1e-3) const;
  // double computeFluxUnm(double E, double zMax, double relError = 1e-3) const;

 protected:
  double dbdE(double E) const;
  double beta(double E) const;
};

}  // namespace solutions
}  // namespace simprop

#endif