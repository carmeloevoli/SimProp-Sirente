#include "simprop.h"

using namespace simprop;

int main(int argc, char* argv[]) {
  try {
    log::startup_information();
    if (argc == 2) {
      utils::Timer timer("SimProp total time");

      Params params(argv[1]);
      params.print();

      SimProp simprop(params);
      simprop.buildInitialStates();
      simprop.dumpPrimaryParticles("output/initial.txt");
      simprop.run();

      // OutputManager output(simprop);
      // output.save();

    } else {
      throw std::runtime_error("please provide an input file as './simprop params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}