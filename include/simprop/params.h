#ifndef SIMPROP_PARAMS_H
#define SIMPROP_PARAMS_H

#include <string>

#include "simprop/Units.h"
#include "simprop/logging.h"
#include "simprop/pid.h"

using EnergyRange = std::pair<double, double>;

namespace simprop {

class Params {
 private:
  std::string m_simName = "fiducial";
  unsigned int m_seed = 12345;
  unsigned int m_nParticles = 1000;
  EnergyRange m_energyRange = {1e9 * SI::eV, 1e20 * SI::eV};
  double m_maxRedshift = 2.0;
  PID m_pid = Fe56;

 public:
  explicit Params(const char* inputFilename);
  virtual ~Params();
  void print();
  const std::string& simName = m_simName;
  const unsigned int& seed = m_seed;
  const unsigned int& nParticles = m_nParticles;
  const EnergyRange& energyRange = m_energyRange;
  const double& maxRedshift = m_maxRedshift;
  const PID& pid = m_pid;

 private:
  //   void set_from_file(const std::string& filename);
  //   void set_params(const std::string& key, const double& value);
};

}  // namespace simprop

#endif  // SIMPROP_PARAMS_H
