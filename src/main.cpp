#include "simprop.h"

using Timer = SimProp::utils::Timings;

int main(int argc, char* argv[]) {
  try {
    SimProp::log::startup_information();
    if (argc == 2) {
      Timer timer;
      timer.StartTiming("Main Program");

      SimProp::Params params;
      params.print();

      //   SOFIA::Sofia sofia(params);

      timer.EndTiming("Main Program");
      timer.PrintAllTimings();
    } else {
      throw std::runtime_error("please provide an input file as './simprop params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}