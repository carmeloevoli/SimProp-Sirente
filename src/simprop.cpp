#include "simprop/simprop.h"

#include <algorithm>

#include "simprop/common.h"

namespace simprop {

SimProp::SimProp(const Params& params)
    : m_params(params), m_rng(RandomNumberGenerator(params.seed)) {
  LOGD << "SimProp constructor";
}

SimProp::~SimProp() { LOGD << "SimProp destructor"; }

void SimProp::buildInitialStates() {
  m_istates.resize(m_params.nParticles);
  std::generate(m_istates.begin(), m_istates.end(), [this]() {
    return std::make_pair(GetRndRedshift(m_params.maxRedshift, 2, this->m_rng()),
                          GetRndEnergy(m_params.energyRange, this->m_rng()));
  });

  // for (auto& p : m_istates) std::cout << p.first << " " << p.second << "\n";
}

void SimProp::run() {}

}  // namespace simprop