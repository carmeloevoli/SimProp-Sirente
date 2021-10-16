#include "simprop/timer.h"

#include "simprop/logging.h"

namespace simprop {
namespace utils {

Timer::Timer() { m_start = AwesomeClock::now(); }

Timer::~Timer() {
  m_end = std::chrono::high_resolution_clock::now();
  m_duration = m_end - m_start;
  LOGI << "SimProp running time : " << m_duration.count() << " s.";
}

}  // namespace utils
}  // namespace simprop
