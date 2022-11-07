#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/core/particle.h"
#include "simprop/core/units.h"

using Range = std::pair<double, double>;

namespace simprop {

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }
inline double energyToWavelenght(double energy) { return SI::hPlanck * SI::cLight / energy; }

double getRndLogUniform(Range range, double r);
Range getRedshiftRange(const ParticleStack& stack);
Range getGammaRange(const ParticleStack& stack);

}  // namespace simprop

#endif