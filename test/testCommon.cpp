#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(Common, energyToFrequency) {
  double photonEnergy = SI::eV;
  double photonFrenquency = energyToFrequency(photonEnergy);
  EXPECT_NEAR(photonFrenquency, 2.418e14, 1e11);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop