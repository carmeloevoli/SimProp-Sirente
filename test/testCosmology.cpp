#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(Cosmology, hubbleTime) {
  auto cosmology = cosmo::Cosmology();
  EXPECT_DOUBLE_EQ(cosmology.OmegaM + cosmology.OmegaL, 1.);
  EXPECT_DOUBLE_EQ(cosmology.hubbleRate(0.), cosmology.H0);
  EXPECT_LT(cosmology.hubbleRate(0.), cosmology.hubbleRate(1.));
  EXPECT_LT(cosmology.hubbleRate(1.), cosmology.hubbleRate(50.));
}

TEST(Cosmology, comovingDistance) {
  auto cosmology = cosmo::Cosmology();
  EXPECT_DOUBLE_EQ(cosmology.redshift2ComovingDistance(0.), 0.);
  EXPECT_LT(cosmology.redshift2ComovingDistance(0.), cosmology.redshift2ComovingDistance(1.));
  EXPECT_LT(cosmology.redshift2ComovingDistance(1.), cosmology.redshift2ComovingDistance(50.));
  EXPECT_THROW(cosmology.redshift2ComovingDistance(-1.), std::invalid_argument);
  EXPECT_THROW(cosmology.redshift2ComovingDistance(101.), std::invalid_argument);
}

TEST(Cosmology, Planck2018) {
  auto cosmology = cosmo::Planck2018();
  EXPECT_DOUBLE_EQ(cosmology.h, 0.674);
  EXPECT_NEAR(1. / cosmology.H0, 14.516915 * SI::Gyr, SI::Myr);
  EXPECT_DOUBLE_EQ(cosmology.OmegaL, 0.685);
  EXPECT_NEAR(cosmology.OmegaB, 0.0492, 1e-4);
  EXPECT_NEAR(cosmology.OmegaC, 0.2642, 1e-4);
  EXPECT_NEAR(cosmology.OmegaM, 0.3134, 1e-4);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop