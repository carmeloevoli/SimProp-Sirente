#include "simprop/SimProp.h"

#include <algorithm>
#include <fstream>

#include "simprop/common.h"
#include "simprop/utils/timer.h"
#include "simprop/utils/tqdm.h"

namespace simprop {

SimProp::SimProp(const Params& params)
    : m_params(params), m_size(params.nParticles), m_rng(RandomNumberGenerator(params.seed)) {
  LOGD << "SimProp constructor";
}

SimProp::~SimProp() { LOGD << "SimProp destructor"; }

void SimProp::printRanges() const {
  auto zrange = std::minmax_element(m_primaries.begin(), m_primaries.end(),
                                    [](const Particle& a, const Particle& b) { return a.z < b.z; });
  LOGD << "redshift range (" << zrange.first->z << "," << zrange.second->z << ")";
  auto Erange = std::minmax_element(m_primaries.begin(), m_primaries.end(),
                                    [](const Particle& a, const Particle& b) { return a.E < b.E; });
  LOGD << "energy range (" << Erange.first->E / SI::GeV << "," << Erange.second->E / SI::GeV << ")";
}

void SimProp::buildInitialStates() {
  m_primaries.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = GetRndRedshift(m_params.maxRedshift, 2, m_rng());
    const auto E_i = GetRndEnergy(m_params.energyRange, m_rng());
    auto p = Particle{m_params.pid, z_i, E_i};
    m_primaries.emplace_back(p);
  }
  LOGD << "built primaries with size " << m_primaries.size();
  printRanges();
}

void SimProp::dumpPrimaryParticles(std::string filename) {
  LOGD << "dumping " << m_primaries.size() << " particles on " << filename;
  std::ofstream ofile;
  ofile.open(filename);
  for (auto it = m_primaries.begin(); it != m_primaries.end(); it++) {
    ofile << *it;
    ofile << std::endl;
  }
  ofile.close();
}

void SimProp::run() {
  int N = 1000;
  tqdm bar;
  LOGI << "Running SimProp : ";
  utils::Timer timer("SimProp core time");
  bar.set_theme_vertical();
  bar.disable_colors();
  for (int i = 0; i < N; i++) {
    bar.progress(i, N);
    usleep(3000);
  }
  bar.finish();
  LOGI << "done!";
}

}  // namespace simprop