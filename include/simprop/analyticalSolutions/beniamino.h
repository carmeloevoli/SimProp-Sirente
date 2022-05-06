#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/cosmology.h"
#include "simprop/energyLosses/AdiabaticContinuousLosses.h"
#include "simprop/energyLosses/BGG2002ContinuousLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

namespace simprop {
namespace solutions {

class Beniamino {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<losses::AdiabaticContinuousLosses> m_ad;
  std::shared_ptr<losses::BGG2002ContinuousLosses> m_pp;
  std::shared_ptr<losses::PhotoPionContinuousLosses> m_pion;

  double m_slope = 2.7;
  double m_sourceEmissivity = SI::GeV / SI::m3 / SI::sec;
  double m_sourceEvolution = 0.;
  double m_sourceMaxRedshift = 6.;
  double m_sourceCutoff = -1.;

 public:
  Beniamino();
  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zMax, double relError = 1e-3) const;
  double dilationFactor(double E, double zMax, double relError = 1e-3) const;
  double computeFlux(double E) const;
  double computeFluxUnm(double E) const;

 public:
  inline void setSlope(const double &slope) { m_slope = slope; };
  inline void setSourceEvolution(const double &m) { m_sourceEvolution = m; };
  inline void setMaxRedshift(const double &z) { m_sourceMaxRedshift = z; };
  inline void setSourceCutoff(const double &E) { m_sourceCutoff = E; };
};

}  // namespace solutions
}  // namespace simprop

#endif