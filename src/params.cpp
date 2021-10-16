#include "simprop/params.h"

//#include <plog/Log.h>

namespace simprop {

Params::Params() { LOGD << "Params constructor"; }

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

void Params::print() {
  LOGI << "model name   : " << simName;
  LOGI << "min Energy   : " << energyRange.first / SI::GeV << " GeV";
  LOGI << "max Energy   : " << energyRange.second / SI::GeV << " GeV";
  LOGI << "max redshift : " << maxRedshift;
}

}  // namespace simprop