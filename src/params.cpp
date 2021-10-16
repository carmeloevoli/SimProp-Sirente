#include "simprop/params.h"

#include <fstream>

#include "simprop/utils.h"
namespace simprop {

Params::Params(const char* inputFilename) {
  LOGD << "Params constructor";
  m_simName = utils::removeExtensionIniFilename(inputFilename);
}

Params::~Params() { LOGD << "Params destructor"; }

// void Params::set_params(const std::string& key, const double& value) {
//   if (key == "H")
//     _H = value * cgs::kpc;
//   else if (key == "dt")
//     _dt = value * cgs::kyr;
//   else if (key == "u")
//     _u = value * cgs::km / cgs::sec;
//   else if (key == "v_A")
//     _v_A = value * cgs::km / cgs::sec;
//   else if (key == "slope")
//     _injection_slope = value;
//   else if (key == "phi")
//     _modulation_potential = value * cgs::GeV;
//   else if (key == "id")
//     _id = (int)value;
// }

// void Params::set_from_file(const std::string& filename) {
//   std::ifstream infile(filename.c_str());
//   std::string line;
//   while (std::getline(infile, line)) {
//     std::istringstream iss(line);
//     std::string key;
//     double value;
//     if (!(iss >> key >> value)) {
//       break;
//     }  // error
//     set_params(key, value);
//   }
// }

#define PFILE(NAME, VALUE) \
  (pfile << std::left << std::setw(20) << std::setfill('_') << NAME << VALUE << std::endl)

void Params::print() {
  std::ofstream pfile;
  pfile.open(simName + ".params");
  PFILE("model name", simName);
  PFILE("seed", seed);
  PFILE("N Particles", nParticles);
  PFILE("min Energy [GeV]", energyRange.first / SI::GeV);
  PFILE("max Energy [GeV]", energyRange.second / SI::GeV);
  PFILE("max redshift []", maxRedshift);
  pfile.close();
}

}  // namespace simprop