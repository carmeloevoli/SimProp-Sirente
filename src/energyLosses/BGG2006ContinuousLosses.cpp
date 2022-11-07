#include "simprop/energyLosses/BGG2006ContinuousLosses.h"

#include "simprop/utils/logging.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace losses {

BGG2006ContinuousLosses::BGG2006ContinuousLosses() : ContinuousLosses() {
  m_totalLosses.loadTable(totalLossesFilename);
  LOGD << "calling " << __func__ << " constructor";
}

double BGG2006ContinuousLosses::getInterpolated(double E) const {
  double b_l = 0;
  const auto logE = std::log10(E / SI::eV);
  if (m_totalLosses.xIsInside(logE)) {
    b_l = std::pow(10., m_totalLosses.get(logE));
    b_l /= SI::year;
  }
  return b_l;
}

double BGG2006ContinuousLosses::beta(PID pid, double Gamma, double z) const {
  // See 10.1016/j.astropartphys.2012.07.010, eq. 3 and 5
  const auto E_over_A = Gamma * SI::protonMassC2;
  const auto redshiftedEnergy = E_over_A * (1. + z);
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
