#ifndef SIMPROP_COSMOLOGY_PLANCK2018_H
#define SIMPROP_COSMOLOGY_PLANCK2018_H

#include "simprop/units.h"

namespace simprop {
namespace planck2018 {

static constexpr double T_CMB = 2.7255 * SI::K;
static constexpr double h = 0.674;
static constexpr double h2 = h * h;
static constexpr double H0 = 100. * h * SI::km / SI::sec / SI::Mpc;
static constexpr double OmegaR = 2.473e-5 / h2;
static constexpr double OmegaB = 0.02237 / h2;
static constexpr double OmegaC = 0.1200 / h2;
static constexpr double OmegaM = OmegaB + OmegaC;
static constexpr double OmegaL = 0.685;

}  // namespace planck2018
}  // namespace simprop

#endif