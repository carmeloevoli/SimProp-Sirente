#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/units.h"

namespace simprop {

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }
inline double energyToWavelenght(double energy) { return SI::hPlanck * SI::cLight / energy; }

double getRndGamma(std::pair<double, double> gammaRange, double slope, double r);
double getRndRedshift(std::pair<double, double> redshiftRange, int evolutionIndex, double r);

}  // namespace simprop

#endif