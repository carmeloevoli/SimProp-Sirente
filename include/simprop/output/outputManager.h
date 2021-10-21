#ifndef SIMPROP_OUTPUT_OUTPUTMANAGER_H
#define SIMPROP_OUTPUT_OUTPUTMANAGER_H

#include "simprop/SimProp.h"

namespace simprop {

class OutputManager {
 public:
  OutputManager(const SimProp& simprop) {}
  virtual ~OutputManager() = default;
  void save() const {};
};

}  // namespace simprop

#endif