// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/energyLosses/AdiabaticContinuousLosses.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace losses {

AdiabaticContinuousLosses::AdiabaticContinuousLosses(
    const std::shared_ptr<cosmo::Cosmology>& cosmology)
    : ContinuousLosses(), m_cosmology(cosmology) {
  LOGD << "calling " << __func__ << " constructor";
}

double AdiabaticContinuousLosses::beta(PID pid, double Gamma, double z) const {
  const auto b_a = 1. / (1. + z);
  return b_a / m_cosmology->dtdz(z);
}

}  // namespace losses
}  // namespace simprop
