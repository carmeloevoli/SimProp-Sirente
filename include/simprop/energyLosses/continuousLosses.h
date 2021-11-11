#ifndef SIMPROP_ENERGYLOSSES_CONTINUOUS_H
#define SIMPROP_ENERGYLOSSES_CONTINUOUS_H

#include "simprop/cosmology/cosmology.h"
#include "simprop/pid.h"
#include "simprop/units.h"
#include "simprop/utils/lookupTable.h"

namespace simprop {
namespace losses {

class ContinuousLosses {
 protected:
  const std::string totalLossesFilename = "data/losses_pair+pion_BGG2002.txt";
  utils::LookupTable<81, 1> m_totalLosses{totalLossesFilename};

 public:
  ContinuousLosses() {}
  virtual ~ContinuousLosses() = default;

  double dlnGamma_dz(double z, double E, PID pid) const;
  double evolve(double E_i, double z_i, double z_f, PID pid) const;
  double evolve_rk4(double E_i, double z_i, double z_f, PID pid) const;
};

}  // namespace losses

}  // namespace simprop

#endif