#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/units.h"

namespace simprop {

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }
inline double energyToWavelenght(double energy) { return SI::hPlanck * SI::cLight / energy; }

}  // namespace simprop

#endif