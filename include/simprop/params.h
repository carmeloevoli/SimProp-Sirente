#ifndef INCLUDE_SIMPROP_PARAMS_H
#define INCLUDE_SIMPROP_PARAMS_H

#include <string>

#include "simprop/Units.h"
#include "simprop/logging.h"

namespace simprop {

using namespace units;

class Params {
 private:
  std::string m_simName = "fiducial";
  QEnergy m_minEnergy = 1e9_eV;
  QEnergy m_maxEnergy = 1e20_eV;
  double m_maxRedshift = 2.0;

 public:
  Params();
  virtual ~Params();
  void print();

  //   void set_from_file(const std::string& filename);
  //   void set_params(const std::string& key, const double& value);

  const std::string& simName = m_simName;
  const QEnergy& minEnergy = m_minEnergy;
  const QEnergy& maxEnergy = m_maxEnergy;
  const double& maxRedshift = m_maxRedshift;
};

}  // namespace simprop

#endif /* INCLUDE_SIMPROP_PARAMS_H */
