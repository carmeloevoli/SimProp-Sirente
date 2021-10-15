#ifndef SIMPROP_UNITSSIMPROP_H
#define SIMPROP_UNITSSIMPROP_H

#include "UnitsDerived.h"

namespace SimProp {
namespace units {

// SIMPROP specific types
// ------------------------------------
//             l  t  m  I  T  N  J  A
QUANTITY_TYPE(2, 0, 0, 1, 0, 0, 0, 0, 0, QSynchroConstant);
QUANTITY_TYPE(-5, 1, -1, 0, 0, 0, 0, 0, 0, QGREmissivity);

}  // namespace units
}  // namespace SimProp

#endif  // SIMPROP_UNITSSIMPROP_H
