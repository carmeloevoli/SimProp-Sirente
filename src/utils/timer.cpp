#include "simprop/utils/timer.h"

#include "simprop/utils/logging.h"

namespace simprop {
namespace utils {

Timer::Timer(std::string message) : m_message(message), m_start(AwesomeClock::now()) {}

Timer::~Timer() {
  m_end = std::chrono::high_resolution_clock::now();
  m_duration = m_end - m_start;
  LOGI << m_message << " " << m_duration.count() << " s.";
}

}  // namespace utils
}  // namespace simprop
