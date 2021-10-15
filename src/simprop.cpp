#include "simprop/simprop.h"

#include <algorithm>
#include <cstdlib>  // std::rand, std::srand

namespace simprop {

double GenerateRandomRedshift() { return (double)(std::rand() % 100); }
QEnergy GenerateRandomEnergy() { return (QEnergy)(std::rand() % 100); }

SimProp::SimProp(const Params& params) : m_rng(utils::RNG<double>(params.seed)) {
  m_redshifts.resize(params.nParticles);
  m_energies.resize(params.nParticles);
}

void SimProp::init() {
  std::generate(m_redshifts.begin(), m_redshifts.end(), GenerateRandomRedshift);
  std::generate(m_energies.begin(), m_energies.end(), GenerateRandomEnergy);
}

void SimProp::run() {}

}  // namespace simprop