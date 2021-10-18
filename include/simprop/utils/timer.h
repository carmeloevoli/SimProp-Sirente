#ifndef SIMPROP_UTILS_TIMER_H
#define SIMPROP_UTILS_TIMER_H

#include <chrono>
#include <string>

namespace simprop {
namespace utils {

using AwesomeClock = std::chrono::high_resolution_clock;

class Timer {
 protected:
  std::string m_message;
  std::chrono::time_point<AwesomeClock> m_start;
  std::chrono::time_point<AwesomeClock> m_end;
  std::chrono::duration<double> m_duration;

 public:
  explicit Timer(std::string message);
  ~Timer();
};

}  // namespace utils
}  // namespace simprop

#endif  // SIMPROP_UTILS_TIMER_H
