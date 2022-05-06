#include "simprop/energyLosses/BGG2002ContinuousLosses.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

BGG2002ContinuousLosses::BGG2002ContinuousLosses() : ContinuousLosses() {
  LOGD << "calling " << __func__ << " constructor";
  m_totalLosses.loadTable(totalLossesFilename);
}

double BGG2002ContinuousLosses::getInterpolated(double E) const {
  double b_l = 0;
  const auto logE = std::log10(E / SI::eV);
  if (m_totalLosses.xIsInside(logE)) {
    b_l = std::pow(10., m_totalLosses.get(logE));
    b_l /= SI::year;
  }
  return b_l;
}

double BGG2002ContinuousLosses::beta(PID pid, double Gamma, double z) const {
  const auto E = Gamma * SI::protonMassC2 * getPidNucleusMassNumber(pid);
  const auto redshiftedEnergy = E * (1. + z);
  double b_l = getInterpolated(redshiftedEnergy);
  if (b_l > 0.) {
    b_l *= pow3(1. + z);
    const double Z = (double)getPidNucleusCharge(pid);
    const double A = (double)getPidNucleusMassNumber(pid);
    b_l *= pow2(Z) / A;
  }
  return std::max(b_l, 0.);
}

}  // namespace losses
}  // namespace simprop
