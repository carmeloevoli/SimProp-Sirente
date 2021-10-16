#include "simprop/simprop.h"

#include <algorithm>

#include "simprop/common.h"

namespace simprop {

SimProp::SimProp(const Params& params)
    : m_params(params), m_size(params.nParticles), m_rng(RandomNumberGenerator(params.seed)) {
  LOGD << "SimProp constructor";
}

SimProp::~SimProp() { LOGD << "SimProp destructor"; }

void SimProp::buildInitialStates() {
  m_primaries.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = GetRndRedshift(m_params.maxRedshift, 2, this->m_rng());
    const auto E_i = GetRndEnergy(m_params.energyRange, this->m_rng());
    auto p = PrimaryParticle(m_params.pid, z_i, E_i);
    m_primaries.emplace_back(p);
  }
}

void SimProp::run() {}

}  // namespace simprop