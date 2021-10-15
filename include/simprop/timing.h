#ifndef INCLUDE_SIMPROP_TIMING_H
#define INCLUDE_SIMPROP_TIMING_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>

namespace SimProp {
namespace utils {

// credits to FML https://fml.wintherscoming.no/

using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

class Timings {
 private:
  std::map<std::string, TimePoint> timings;
  std::map<std::string, double> elapsed_time_sec;
  std::mutex timings_mutex;

  TimePoint getTime() { return std::chrono::steady_clock::now(); }

  // Convert an interval to seconds
  double timeInSeconds(const TimePoint& time_start, const TimePoint& time_end) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(time_end - time_start).count();
  }

 public:
  Timings() = default;

  /// Output all the recorded timings
  void PrintAllTimings() {
    for (const auto& t : elapsed_time_sec) {
      LOGI << "Total elapsed time for [" << std::setw(35) << t.first << "]: " << std::setw(10) << t.second << " sec";
    }
  }

  /// Start timing
  /// @param[in] name The label to give to the timing
  ///
  void StartTiming(const char* name) { StartTiming(std::string(name)); }

  /// Start timing
  void StartTiming(std::string name) {
    auto start_time = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> guard(timings_mutex);

    auto it = timings.find(name);
    if (it != timings.end()) {
      it->second = start_time;
    } else {
      timings[name] = start_time;
      elapsed_time_sec[name] = 0.0;
    }
  }

  /// End timing
  /// @param[in] name The label to give to the timing
  /// @param[in] print Optional: print timing to screen right away
  ///
  double EndTiming(const char* name, bool print = false) { return EndTiming(std::string(name), print); }

  /// End timing and print to screen if print = true
  /// @param[in] name The label to give to the timing
  /// @param[in] print Optional: print timing to screen right away
  ///
  double EndTiming(std::string name, bool print = false) {
    auto end_time = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> guard(timings_mutex);

    double time_sec = 0.0;
    double tot_time_sec = 0.0;

    // Find current time
    auto it = timings.find(name);
    if (it != timings.end()) {
      auto start_time = it->second;
      time_sec = timeInSeconds(start_time, end_time);
    }

    // Find total elapsed time
    auto it2 = elapsed_time_sec.find(name);
    if (it2 == elapsed_time_sec.end()) {
      elapsed_time_sec[name] = time_sec;
    } else {
      it2->second += time_sec;
      tot_time_sec = it2->second;
    }

    if (print) {
      if (tot_time_sec == 0.0) {
        if (time_sec == 0.0) {
          LOGE << "Can't print elapsed time. Start time for time-point [" << name << "] was not found";
        } else {
          LOGI << "Elapsed time for [" << name << "]: " << time_sec << " sec";
        }
      } else {
        LOGI << "Current elapsed time for [" << name << "]: " << time_sec << " sec. Total: " << tot_time_sec << " sec";
      }
    }

    return time_sec;
  }

  /// Print to screen total time for a given label
  /// @param[in] name The label to print the time of
  ///
  void PrintTotalTime(const char* name) { PrintTotalTime(std::string(name)); }

  /// Print to screen total time for a given label
  /// @param[in] name The label to print the time of
  ///
  void PrintTotalTime(std::string name) {
    auto it = elapsed_time_sec.find(name);
    if (it != elapsed_time_sec.end()) {
      double time_sec = it->second;
      LOGI << "Total elapsed time for [" << name << "]: " << time_sec << " sec";
    }
  }
};

}  // namespace utils
}  // namespace SimProp

#endif