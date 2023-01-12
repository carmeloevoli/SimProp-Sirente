#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/core/cosmology.h"
#include "simprop/energyLosses/BGG2006ContinuousLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

namespace simprop {
namespace solutions {

class Beniamino {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<losses::BGG2006ContinuousLosses> m_pp;
  std::shared_ptr<losses::PhotoPionContinuousLosses> m_pion;

  const double m_sourceEmissivity = SI::GeV / SI::m3 / SI::sec;
  const double m_maxEnergy = 1e23 * SI::eV;
  double m_slope = 2.7;
  double m_sourceEvolution = 0.;
  double m_sourceCutoff = -1.;
  bool m_doPhotoPion = false;

 public:
  Beniamino();
  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zMax, double relError = 1e-3) const;
  double dilationFactor(double E, double zMax, double relError = 1e-3) const;
  double computeFlux(double E, double zMax, double relError = 1e-3) const;
  double computeFluxUnm(double E, double zMax, double relError = 1e-3) const;

 public:
  inline void setSlope(const double &slope) { m_slope = slope; };
  inline void setSourceEvolution(const double &m) { m_sourceEvolution = m; };
  inline void setSourceCutoff(const double &E) { m_sourceCutoff = E; };
  inline void disablePhotoPion() { m_doPhotoPion = false; }
  inline void enablePhotoPion() { m_doPhotoPion = true; }

 protected:
  double dbdE(double E) const;
};

}  // namespace solutions
}  // namespace simprop

#endif