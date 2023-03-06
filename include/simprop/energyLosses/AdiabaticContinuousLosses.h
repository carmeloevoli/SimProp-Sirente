// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_ENERGYLOSSES_ADIABATICCONTINUOUSLOSSES_H_
#define SIMPROP_ENERGYLOSSES_ADIABATICCONTINUOUSLOSSES_H_

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"

namespace simprop {
namespace losses {

class AdiabaticContinuousLosses final : public ContinuousLosses {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;

 public:
  AdiabaticContinuousLosses(const std::shared_ptr<cosmo::Cosmology>& cosmology);
  virtual ~AdiabaticContinuousLosses() = default;

  double beta(PID pid, double Gamma, double z = 0) const override;
};

}  // namespace losses
}  // namespace simprop

#endif  // SIMPROP_ENERGYLOSSES_ADIABATICCONTINUOUSLOSSES_H_
