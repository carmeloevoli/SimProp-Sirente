#ifndef INCLUDE_SIMPROP_PARAMS_H
#define INCLUDE_SIMPROP_PARAMS_H

#include <string>

#include "include/simprop.h"

//#include "simprop/Units.h"

namespace SimProp {

class Params {
 private:
  std::string m_simName = "fiducial";
  // QEnergy minEnergy = 1e9_GeV;
  // QEnergy maxEnergy = 1e20_GeV;

 public:
  Params();
  virtual ~Params();
  void print();

  //   void set_from_file(const std::string& filename);
  //   void set_params(const std::string& key, const double& value);

  const std::string& simName = m_simName;
  // const Wavenumber& kMin = m_kMin;
  // const Wavenumber& kMax = m_kMax;
  // const Length& H = m_haloSize;
};

}  // namespace SimProp

#endif /* INCLUDE_SIMPROP_PARAMS_H */
