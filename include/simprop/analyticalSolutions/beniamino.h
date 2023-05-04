#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/core/cosmology.h"
#include "simprop/energyLosses/PairProductionLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

namespace simprop {
namespace solutions {

class Beniamino {
 public:
  Beniamino(bool doPhotoPion = false);
  Beniamino(double injSlope, double sourceEvolution, double sourceCutoff, bool doPhotoPion = false);
  Beniamino& doCaching();

  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zNow, double zMax, double relError = 1e-3) const;
  double dilationFactor(double E, double zNow, double zMax, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double zMax, double relError = 1e-2) const;
  // double computeFluxUnm(double E, double zMax, double relError = 1e-3) const;

 public:
  double dbdE(double E) const;
  double beta(double E) const;

 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<losses::PairProductionLosses> m_pair;
  std::shared_ptr<losses::PhotoPionContinuousLosses> m_pion;
  utils::LookupArray<8000> m_losses;

  const double m_sourceEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  const double m_maxEnergy{1e23 * SI::eV};
  const double m_minEnergy{1e17 * SI::eV};

  double m_slope{2.6};
  double m_sourceEvolution{0.};
  double m_sourceCutoff{-1.};
  bool m_doPhotoPion{false};
  bool m_doCaching{false};
};

}  // namespace solutions
}  // namespace simprop

#endif