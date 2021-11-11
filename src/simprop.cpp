#include "simprop/SimProp.h"

#include <algorithm>
#include <fstream>

#include "simprop/common.h"
#include "simprop/cosmology/cosmology.h"
#include "simprop/energyLosses/BGG2002ContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/Dominguez2011PhotonField.h"
#include "simprop/photonFields/Gilmore2012PhotonField.h"
#include "simprop/utils/misc.h"
#include "simprop/utils/timer.h"
#include "simprop/utils/tqdm.h"

namespace simprop {

SimProp::SimProp(const Params& params)
    : m_params(params), m_size(params.nParticles), m_rng(RandomNumberGenerator(params.seed)) {
  LOGD << "SimProp constructor";
}

SimProp::~SimProp() { LOGD << "SimProp destructor"; }

void SimProp::buildInitialStates() {
  m_particles.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = GetRndRedshift(m_params.redshiftRange, 2, m_rng());
    const auto E_i = GetRndEnergy(m_params.energyRange, 1, m_rng());
    auto p = Particle{m_params.pid, z_i, E_i};
    m_particles.emplace_back(p);
  }
  LOGD << "built primaries with size " << m_particles.size();
  printStateRanges();
}

void SimProp::printStateRanges() const {
  auto z_r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getRedshift() < b.getRedshift(); });
  LOGD << "z range (" << z_r.first->getRedshift() << "," << z_r.second->getRedshift() << ")";
  auto E_r = std::minmax_element(
      m_particles.begin(), m_particles.end(),
      [](const Particle& a, const Particle& b) { return a.getEnergy() < b.getEnergy(); });
  using SI::eV;
  LOGD << "E range (" << E_r.first->getEnergy() / eV << "," << E_r.second->getEnergy() / eV << ")";
}

void SimProp::dumpParticles(std::string filename) const {
  LOGD << "dumping " << m_particles.size() << " particles on " << filename;
  std::ofstream ofile;
  ofile.open(filename);
  for (auto& p : m_particles) {
    ofile << p;
    ofile << std::endl;
  }
  ofile.close();
}

void SimProp::buildPhotonFields() {
  m_photonFields.emplace_back(new photonfield::CMB());
  // if (m_params.eblModel == Params::GILMORE2012)
  //   m_photonFields.emplace_back(new photonfield::Gilmore2012Field());
  // else if (m_params.eblModel == Params::DOMINGUEZ2011)
  //   m_photonFields.emplace_back(new photonfield::Dominguez2011Field());
}

void SimProp::buildContinuousLosses() {
  m_continuousLosses = std::make_shared<losses::BGG2002ContinuousLosses>();
}

bool isActive(Particle p) {
  const double minPropagatingEnergy = SI::eV;
  return (p.getRedshift() > 1e-20 && p.getEnergy() > minPropagatingEnergy);
}

void SimProp::Evolve(ParticleStack::iterator it) {
  const double dz = 0.01;
  const auto nowRedshift = it->getNow().z;
  const auto nextRedshift = std::max(nowRedshift - dz, 0.);
  it->getNow().z = nextRedshift;
  const auto nowEnergy = it->getNow().E;
  it->getNow().E = m_continuousLosses->evolve(nowEnergy, nowRedshift, nextRedshift, it->getPid());
}

void SimProp::run() {
  LOGI << "Running SimProp : ";
  utils::Timer timer("SimProp core time");
  ParticleStack::iterator it = m_particles.begin();
  while (it != m_particles.end()) {
    Evolve(it);
    it = std::find_if(m_particles.begin(), m_particles.end(), isActive);
  }
  LOGI << "done!";
}

}  // namespace simprop