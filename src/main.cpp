#include "simprop.h"

using namespace SimProp;

int main(int argc, char* argv[]) {
  try {
    SimProp::log::startup_information();
    if (argc == 2) {
      utils::Timer timer;

      Params params;
      params.print();

      //   SOFIA::Sofia sofia(params);
    } else {
      throw std::runtime_error("please provide an input file as './simprop params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}