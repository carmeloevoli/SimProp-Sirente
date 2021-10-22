#ifndef SIMPROP_PARAMS_H
#define SIMPROP_PARAMS_H

#include <string>

#include "simprop/Units.h"
#include "simprop/pid.h"
#include "simprop/utils/logging.h"

using Range = std::pair<double, double>;

namespace simprop {

class Params {
 public:
  enum EblModel { GILMORE2012, DOMINGUEZ2011 };
  std::string toString(EblModel model);

 private:
  std::string m_simName = "fiducial";
  unsigned int m_seed = 12345;
  unsigned int m_nParticles = 1000;
  Range m_energyRange = {1e9 * SI::eV, 1e20 * SI::eV};
  Range m_redshiftRange = {0., 2.};
  PID m_pid = Fe56;
  EblModel m_eblModel = GILMORE2012;

 public:
  explicit Params(const char* inputFilename);
  virtual ~Params();
  void print();
  const std::string& simName = m_simName;
  const unsigned int& seed = m_seed;
  const unsigned int& nParticles = m_nParticles;
  const Range& energyRange = m_energyRange;
  const Range& redshiftRange = m_redshiftRange;
  const PID& pid = m_pid;
  const EblModel& eblModel = m_eblModel;

 private:
  //   void set_from_file(const std::string& filename);
  //   void set_params(const std::string& key, const double& value);
};

}  // namespace simprop

#endif  // SIMPROP_PARAMS_H
