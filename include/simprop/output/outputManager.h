#ifndef SIMPROP_OUTPUT_OUTPUTMANAGER_H
#define SIMPROP_OUTPUT_OUTPUTMANAGER_H

#include "simprop/simprop.h"

namespace simprop {

class OutputManager {
 public:
  OutputManager(const SimProp& simprop) {}
  virtual ~OutputManager() = default;
  void save() const {};
};

}  // namespace simprop

#endif