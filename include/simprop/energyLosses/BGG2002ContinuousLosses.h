#ifndef SIMPROP_LOSSES_BGG2002_CONTINUOUS_H
#define SIMPROP_LOSSES_BGG2002_CONTINUOUS_H

#include "simprop/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/units.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace losses {

class BGG2002ContinuousLosses final : public ContinuousLosses {
 protected:
  const std::string totalLossesFilename = "data/losses_pair_BGG2002.txt";
  utils::LookupArray<501> m_totalLosses{totalLossesFilename};

 public:
  explicit BGG2002ContinuousLosses(const cosmo::Cosmology& cosmology)
      : ContinuousLosses(cosmology) {}
  virtual ~BGG2002ContinuousLosses() = default;

  double dlnE_dt(PID pid, double E, double z = 0) const override;
  double dlnE_dz(PID pid, double E, double z = 0) const override;
  // double evolve(double E_i, double z_i, double z_f, PID pid) const override;
  // double evolve_rk4(double E_i, double z_i, double z_f, PID pid) const;

 protected:
  double getInterpolated(double E) const;
};

}  // namespace losses

}  // namespace simprop

#endif