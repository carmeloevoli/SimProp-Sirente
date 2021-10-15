#ifndef INCLUDE_SIMPROP_TIMING_H
#define INCLUDE_SIMPROP_TIMING_H

#include <chrono>

namespace SimProp {
namespace utils {

using AwesomeClock = std::chrono::high_resolution_clock;

class Timer {
 protected:
  std::chrono::time_point<AwesomeClock> m_start;
  std::chrono::time_point<AwesomeClock> m_end;
  std::chrono::duration<double> m_duration;

 public:
  Timer();
  ~Timer();
};

}  // namespace utils
}  // namespace SimProp

#endif  // INCLUDE_SIMPROP_TIMING_H
