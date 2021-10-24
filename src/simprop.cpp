#include "simprop/SimProp.h"

#include <algorithm>
#include <fstream>

#include "simprop/common.h"
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

void SimProp::printRanges() const {
  auto zrange = std::minmax_element(m_primaries.begin(), m_primaries.end(),
                                    [](const Particle& a, const Particle& b) { return a.z < b.z; });
  LOGD << "z initial range (" << zrange.first->z << "," << zrange.second->z << ")";
  auto Erange = std::minmax_element(m_primaries.begin(), m_primaries.end(),
                                    [](const Particle& a, const Particle& b) { return a.E < b.E; });
  using SI::GeV;
  LOGD << "E initial range (" << Erange.first->E / GeV << "," << Erange.second->E / GeV << ")";
}

void SimProp::buildInitialStates() {
  m_primaries.reserve(m_size);
  for (size_t i = 0; i < m_size; ++i) {
    const auto z_i = GetRndRedshift(m_params.redshiftRange.second, 2, m_rng());
    const auto E_i = GetRndEnergy(m_params.energyRange, m_rng());
    auto p = Particle{m_params.pid, z_i, E_i};
    m_primaries.emplace_back(p);
  }
  LOGD << "built primaries with size " << m_primaries.size();
  printRanges();
}

void SimProp::dumpPrimaryParticles() {
  std::string filename = "output/" + m_params.simName + ".primaries";
  LOGD << "dumping " << m_primaries.size() << " particles on " << filename;
  std::ofstream ofile;
  ofile.open(filename);
  for (auto& p : m_primaries) {
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

void SimProp::dumpPhotonFields() {
  std::string filename = "output/" + m_params.simName + ".ebl";
  LOGD << "dumping photon fields on " << filename;
  const auto ePhoton = utils::LogAxis(1e-6 * SI::eV, 1e6 * SI::eV, 1200);
  const auto units = SI::nW / utils::pow<2>(SI::meter) / SI::sr;
  std::ofstream ofile;
  ofile.open(filename);
  ofile << "# energy [eV] - frequency [GHz] - brightness [nW/m^2/sr]\n";
  for (auto e : ePhoton) {
    ofile << e / SI::eV << " ";
    ofile << energyToFrequency(e) / SI::GHz << " ";
    auto n_gamma = std::accumulate(
        m_photonFields.begin(), m_photonFields.end(), 0.,
        [e](double result, const auto& field) { return result + field->getPhotonDensity(e); });
    ofile << utils::pow<2>(e) * n_gamma * SI::cLight / (4. * M_PI) / units;
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