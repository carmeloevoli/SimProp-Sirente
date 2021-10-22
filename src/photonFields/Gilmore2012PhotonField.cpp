#include "simprop/photonFields/Gilmore2012PhotonField.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "simprop/utils/misc.h"

namespace simprop {
namespace photonfield {

Gilmore2012Field::Gilmore2012Field() {
  //   m_fieldName = "Gilmore2012";
  //   if (!utils::fileExist(m_datafile))
  //     throw std::runtime_error("datafile for Gilmore2012 field not found");
  //   std::ifstream datafile(m_datafile.c_str());
  //   const size_t nHeaderLines = 2;
  //   size_t count = 0;
  //   std::string line;
  //   while (getline(datafile, line)) {
  //     if (count >= nHeaderLines) {
  //       // std::cout << line << "\n";
  //       auto items = utils::split(line, " ");
  //       m_ePhoton.push_back(std::log10(atof(items[0].c_str())));
  //       //   //   auto s = items.at(useCol);
  //       //   v.push_back(atof(s.c_str()));
  //       //       d = np.genfromtxt(self.files, unpack=True)
  //       // eps = h * c0 / (d[0] * 1e-10)  # [J]
  //       // n = d[1:] * erg / 1e-4 * d[0] / eps**2 * (4 * np.pi / c0)
  //       // for i,z in enumerate(self.redshift):
  //       //     n[i] /= (1 + z)**3  # make comoving
  //       //     self.data[z] = eps[::-1], n[i][::-1]
  //     }
  //     count++;
  //   }
  //   if (m_ePhoton.size() != 101) throw std::runtime_error("problem in reading file");

  //   datafile.close();
}

double Gilmore2012Field::getPhotonDensity(double ePhoton, double z) const { return 0; };

}  // namespace photonfield
}  // namespace simprop