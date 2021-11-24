#ifndef SIMPROP_INTERACTIONS_ABSTRACTPHOTODISINTEGRATION_H
#define SIMPROP_INTERACTIONS_ABSTRACTPHOTODISINTEGRATION_H

#include "simprop/pid.h"

namespace simprop {
namespace interactions {

class AbstractPhotoDisintegration {
 public:
  AbstractPhotoDisintegration() {}
  virtual ~AbstractPhotoDisintegration() = default;

  virtual double getAbsorptionSigma(double photonEnergy) const = 0;
};

}  // namespace interactions
}  // namespace simprop

#endif