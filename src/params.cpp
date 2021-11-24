#include "simprop/params.h"

#include <fstream>

#include "simprop/utils/misc.h"

namespace simprop {

Params::Params(const char* inputFilename) {
  LOGD << "Params constructor";
  m_simName = utils::removeExtensionIniFilename(inputFilename);
}

Params::~Params() { LOGD << "Params destructor"; }

std::string Params::toString(EblModel model) {
  if (model == GILMORE2012)
    return "Gilmore2012";
  else if (model == DOMINGUEZ2011)
    return "Dominguez2011";
  else
    return "CmbOnly";
}

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

template <typename T>
std::string printElement(T t, const std::string& name) {
  std::stringstream ss;
  ss << std::left << std::setw(20) << name;
  ss << std::left << std::setw(10) << t;
  ss << "\n";
  return ss.str();
}

void Params::print() {
  std::ofstream pfile;
  pfile.open(simName + ".params");
  pfile << printElement(simName, "Model name");
  pfile << printElement(seed, "Seed");
  pfile << printElement(getPidName(pid), "Mass");
  pfile << printElement(nParticles, "Size Particles");
  pfile << printElement(energyRange.first / SI::GeV, "Min Energy [GeV]");
  pfile << printElement(energyRange.second / SI::GeV, "Max Energy [GeV]");
  pfile << printElement(redshiftRange.second, "Max redshift");
  pfile << printElement(toString(eblModel), "EBL model");
  pfile.close();
}

}  // namespace simprop