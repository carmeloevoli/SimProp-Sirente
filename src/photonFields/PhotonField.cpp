// Copyright 2023 SimProp-dev [MIT License]
#include "simprop/photonFields/PhotonField.h"

#include <cmath>

#include "simprop/utils/numeric.h"

namespace simprop {
namespace photonfields {

double PhotonField::computeIntegratedDensity(double z) const {
  auto lnEpsMin = std::log(getMinPhotonEnergy());
  auto lnEpsMax = std::log(getMaxPhotonEnergy());
  auto value = utils::simpsonIntegration<double>(
      [this, z](double lnEps) {
        auto eps = std::exp(lnEps);
        return eps * this->density(eps, z);
      },
      lnEpsMin, lnEpsMax, 1000);
  return value;
}

}  // namespace photonfields
}  // namespace simprop
