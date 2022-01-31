#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(RNG, rangeUnity) {
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  EXPECT_DOUBLE_EQ(rng.min(), 0.);
  EXPECT_DOUBLE_EQ(rng.max(), 1.);
}

TEST(RNG, withinRange) {
  RandomNumberGenerator rng = utils::RNG<double>(5678);
  for (size_t i = 0; i < 1000000; ++i) {
    auto r = rng();
    EXPECT_GE(r, 0.0);
    EXPECT_LE(r, 1.0);
  }
}

TEST(RNG, uniform) {
  RandomNumberGenerator rng = utils::RNG<double>(5678);
  for (size_t i = 0; i < 1000000; ++i) {
    auto r = rng.uniform(-0.1, 0.1);
    EXPECT_GE(r, -0.1);
    EXPECT_LE(r, 0.1);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop