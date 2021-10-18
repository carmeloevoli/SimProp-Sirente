#ifndef SIMPROP_UTILS_TIMER_H
#define SIMPROP_UTILS_TIMER_H

#include <chrono>

namespace simprop {
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
}  // namespace simprop

#endif  // SIMPROP_UTILS_TIMER_H
