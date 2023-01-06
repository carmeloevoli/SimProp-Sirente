#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(Pid, isNucleon) {
  EXPECT_TRUE(pidIsNucleon(proton));
  EXPECT_TRUE(pidIsNucleon(neutron));
  EXPECT_FALSE(pidIsNucleon(He4));
  EXPECT_FALSE(pidIsNucleon(N14));
  EXPECT_TRUE(pidIsNucleus(proton));
  EXPECT_TRUE(pidIsNucleus(neutron));
  EXPECT_TRUE(pidIsNucleus(He4));
  EXPECT_TRUE(pidIsNucleus(N14));
  EXPECT_TRUE(C12 == C12);
  EXPECT_FALSE(neutron == proton);
}

TEST(Pid, removeNucleon) {
  EXPECT_TRUE(removeNucleon(B11, neutron) == B10);
  EXPECT_TRUE(removeNucleon(C12, proton) == B11);
  EXPECT_TRUE(removeNucleon(Ne22, neutron) == Ne21);
  EXPECT_TRUE(removeNucleon(Ti46, proton) == Sc45);
}

TEST(Pid, elementaryName) {
  EXPECT_STREQ("photon", getPidName(photon).c_str());
  EXPECT_STREQ("nu_e", getPidName(neutrino_e).c_str());
  EXPECT_STREQ("antinu_e", getPidName(antineutrino_e).c_str());
  EXPECT_STREQ("nu_mu", getPidName(neutrino_mu).c_str());
  EXPECT_STREQ("antinu_mu", getPidName(antineutrino_mu).c_str());
  EXPECT_STREQ("electron", getPidName(electron).c_str());
  EXPECT_STREQ("positron", getPidName(positron).c_str());
  EXPECT_STREQ("pion_0", getPidName(pionNeutral).c_str());
  EXPECT_STREQ("pion_plus", getPidName(pionPlus).c_str());
  EXPECT_STREQ("pion_minus", getPidName(pionMinus).c_str());
}

TEST(Pid, elementaryNotNucleus) {
  EXPECT_FALSE(pidIsNucleus(photon));
  EXPECT_FALSE(pidIsNucleus(neutrino_e));
  EXPECT_FALSE(pidIsNucleus(antineutrino_e));
  EXPECT_FALSE(pidIsNucleus(neutrino_mu));
  EXPECT_FALSE(pidIsNucleus(antineutrino_mu));
  EXPECT_FALSE(pidIsNucleus(electron));
  EXPECT_FALSE(pidIsNucleus(positron));
  EXPECT_FALSE(pidIsNucleus(pionNeutral));
  EXPECT_FALSE(pidIsNucleus(pionPlus));
  EXPECT_FALSE(pidIsNucleus(pionMinus));
  EXPECT_TRUE(pidIsNucleus(proton));
  EXPECT_TRUE(pidIsNucleus(neutron));
  EXPECT_TRUE(pidIsNucleus(antiproton));
  EXPECT_TRUE(pidIsNucleus(deuterium));
}

TEST(Pid, nucleiName) {
  EXPECT_STREQ("proton", getPidName(proton).c_str());
  EXPECT_STREQ("neutron", getPidName(neutron).c_str());
  EXPECT_STREQ("antiproton", getPidName(antiproton).c_str());
  EXPECT_STREQ("deuterium", getPidName(deuterium).c_str());
  EXPECT_STREQ("He3", getPidName(He3).c_str());
  EXPECT_STREQ("He4", getPidName(He4).c_str());
  EXPECT_STREQ("Be9", getPidName(Be9).c_str());
  EXPECT_STREQ("B10", getPidName(B10).c_str());
  EXPECT_STREQ("B11", getPidName(B11).c_str());
  EXPECT_STREQ("C12", getPidName(C12).c_str());
  EXPECT_STREQ("C13", getPidName(C13).c_str());
  EXPECT_STREQ("N14", getPidName(N14).c_str());
  EXPECT_STREQ("N15", getPidName(N15).c_str());
  EXPECT_STREQ("O16", getPidName(O16).c_str());
  EXPECT_STREQ("O17", getPidName(O17).c_str());
  EXPECT_STREQ("O18", getPidName(O18).c_str());
  EXPECT_STREQ("F19", getPidName(F19).c_str());
  EXPECT_STREQ("Ne20", getPidName(Ne20).c_str());
  EXPECT_STREQ("Ne21", getPidName(Ne21).c_str());
  EXPECT_STREQ("Ne22", getPidName(Ne22).c_str());
  EXPECT_STREQ("Na23", getPidName(Na23).c_str());
  EXPECT_STREQ("Mg24", getPidName(Mg24).c_str());
  EXPECT_STREQ("Mg25", getPidName(Mg25).c_str());
  EXPECT_STREQ("Mg26", getPidName(Mg26).c_str());
  EXPECT_STREQ("Al27", getPidName(Al27).c_str());
  EXPECT_STREQ("Si28", getPidName(Si28).c_str());
  EXPECT_STREQ("Si29", getPidName(Si29).c_str());
  EXPECT_STREQ("Si30", getPidName(Si30).c_str());
  EXPECT_STREQ("P31", getPidName(P31).c_str());
  EXPECT_STREQ("S32", getPidName(S32).c_str());
  EXPECT_STREQ("S33", getPidName(S33).c_str());
  EXPECT_STREQ("S34", getPidName(S34).c_str());
  EXPECT_STREQ("Cl35", getPidName(Cl35).c_str());
  EXPECT_STREQ("Ar36", getPidName(Ar36).c_str());
  EXPECT_STREQ("Cl37", getPidName(Cl37).c_str());
  EXPECT_STREQ("Ar38", getPidName(Ar38).c_str());
  EXPECT_STREQ("K39", getPidName(K39).c_str());
  EXPECT_STREQ("Ca40", getPidName(Ca40).c_str());
  EXPECT_STREQ("K41", getPidName(K41).c_str());
  EXPECT_STREQ("Ca42", getPidName(Ca42).c_str());
  EXPECT_STREQ("Ca43", getPidName(Ca43).c_str());
  EXPECT_STREQ("Ca44", getPidName(Ca44).c_str());
  EXPECT_STREQ("Sc45", getPidName(Sc45).c_str());
  EXPECT_STREQ("Ti46", getPidName(Ti46).c_str());
  EXPECT_STREQ("Ti47", getPidName(Ti47).c_str());
  EXPECT_STREQ("Ti48", getPidName(Ti48).c_str());
  EXPECT_STREQ("Ti49", getPidName(Ti49).c_str());
  EXPECT_STREQ("Cr50", getPidName(Cr50).c_str());
  EXPECT_STREQ("V51", getPidName(V51).c_str());
  EXPECT_STREQ("Cr52", getPidName(Cr52).c_str());
  EXPECT_STREQ("Cr53", getPidName(Cr53).c_str());
  EXPECT_STREQ("Fe54", getPidName(Fe54).c_str());
  EXPECT_STREQ("Mn55", getPidName(Mn55).c_str());
  EXPECT_STREQ("Fe56", getPidName(Fe56).c_str());
}

TEST(Pid, nucleiMass) {
  EXPECT_EQ(1, getPidNucleusMassNumber(proton));
  EXPECT_EQ(1, getPidNucleusMassNumber(neutron));
  EXPECT_EQ(1, getPidNucleusMassNumber(antiproton));
  EXPECT_EQ(2, getPidNucleusMassNumber(deuterium));
  EXPECT_EQ(3, getPidNucleusMassNumber(He3));
  EXPECT_EQ(4, getPidNucleusMassNumber(He4));
  EXPECT_EQ(9, getPidNucleusMassNumber(Be9));
  EXPECT_EQ(10, getPidNucleusMassNumber(B10));
  EXPECT_EQ(11, getPidNucleusMassNumber(B11));
  EXPECT_EQ(12, getPidNucleusMassNumber(C12));
  EXPECT_EQ(13, getPidNucleusMassNumber(C13));
  EXPECT_EQ(14, getPidNucleusMassNumber(N14));
  EXPECT_EQ(15, getPidNucleusMassNumber(N15));
  EXPECT_EQ(16, getPidNucleusMassNumber(O16));
  EXPECT_EQ(17, getPidNucleusMassNumber(O17));
  EXPECT_EQ(18, getPidNucleusMassNumber(O18));
  EXPECT_EQ(19, getPidNucleusMassNumber(F19));
  EXPECT_EQ(20, getPidNucleusMassNumber(Ne20));
  EXPECT_EQ(21, getPidNucleusMassNumber(Ne21));
  EXPECT_EQ(22, getPidNucleusMassNumber(Ne22));
  EXPECT_EQ(23, getPidNucleusMassNumber(Na23));
  EXPECT_EQ(24, getPidNucleusMassNumber(Mg24));
  EXPECT_EQ(25, getPidNucleusMassNumber(Mg25));
  EXPECT_EQ(26, getPidNucleusMassNumber(Mg26));
  EXPECT_EQ(27, getPidNucleusMassNumber(Al27));
  EXPECT_EQ(28, getPidNucleusMassNumber(Si28));
  EXPECT_EQ(29, getPidNucleusMassNumber(Si29));
  EXPECT_EQ(30, getPidNucleusMassNumber(Si30));
  EXPECT_EQ(31, getPidNucleusMassNumber(P31));
  EXPECT_EQ(32, getPidNucleusMassNumber(S32));
  EXPECT_EQ(33, getPidNucleusMassNumber(S33));
  EXPECT_EQ(34, getPidNucleusMassNumber(S34));
  EXPECT_EQ(35, getPidNucleusMassNumber(Cl35));
  EXPECT_EQ(36, getPidNucleusMassNumber(Ar36));
  EXPECT_EQ(37, getPidNucleusMassNumber(Cl37));
  EXPECT_EQ(38, getPidNucleusMassNumber(Ar38));
  EXPECT_EQ(39, getPidNucleusMassNumber(K39));
  EXPECT_EQ(40, getPidNucleusMassNumber(Ca40));
  EXPECT_EQ(41, getPidNucleusMassNumber(K41));
  EXPECT_EQ(42, getPidNucleusMassNumber(Ca42));
  EXPECT_EQ(43, getPidNucleusMassNumber(Ca43));
  EXPECT_EQ(44, getPidNucleusMassNumber(Ca44));
  EXPECT_EQ(45, getPidNucleusMassNumber(Sc45));
  EXPECT_EQ(46, getPidNucleusMassNumber(Ti46));
  EXPECT_EQ(47, getPidNucleusMassNumber(Ti47));
  EXPECT_EQ(48, getPidNucleusMassNumber(Ti48));
  EXPECT_EQ(49, getPidNucleusMassNumber(Ti49));
  EXPECT_EQ(50, getPidNucleusMassNumber(Cr50));
  EXPECT_EQ(51, getPidNucleusMassNumber(V51));
  EXPECT_EQ(52, getPidNucleusMassNumber(Cr52));
  EXPECT_EQ(53, getPidNucleusMassNumber(Cr53));
  EXPECT_EQ(54, getPidNucleusMassNumber(Fe54));
  EXPECT_EQ(55, getPidNucleusMassNumber(Mn55));
  EXPECT_EQ(56, getPidNucleusMassNumber(Fe56));
}

TEST(Pid, nucleiCharge) {
  EXPECT_EQ(1, getPidNucleusCharge(proton));
  EXPECT_EQ(0, getPidNucleusCharge(neutron));
  EXPECT_EQ(-1, getPidNucleusCharge(antiproton));
  EXPECT_EQ(1, getPidNucleusCharge(deuterium));
  EXPECT_EQ(2, getPidNucleusCharge(He3));
  EXPECT_EQ(2, getPidNucleusCharge(He4));
  EXPECT_EQ(4, getPidNucleusCharge(Be9));
  EXPECT_EQ(5, getPidNucleusCharge(B10));
  EXPECT_EQ(5, getPidNucleusCharge(B11));
  EXPECT_EQ(6, getPidNucleusCharge(C12));
  EXPECT_EQ(6, getPidNucleusCharge(C13));
  EXPECT_EQ(7, getPidNucleusCharge(N14));
  EXPECT_EQ(7, getPidNucleusCharge(N15));
  EXPECT_EQ(8, getPidNucleusCharge(O16));
  EXPECT_EQ(8, getPidNucleusCharge(O17));
  EXPECT_EQ(8, getPidNucleusCharge(O18));
  EXPECT_EQ(9, getPidNucleusCharge(F19));
  EXPECT_EQ(10, getPidNucleusCharge(Ne20));
  EXPECT_EQ(10, getPidNucleusCharge(Ne21));
  EXPECT_EQ(10, getPidNucleusCharge(Ne22));
  EXPECT_EQ(11, getPidNucleusCharge(Na23));
  EXPECT_EQ(12, getPidNucleusCharge(Mg24));
  EXPECT_EQ(12, getPidNucleusCharge(Mg25));
  EXPECT_EQ(12, getPidNucleusCharge(Mg26));
  EXPECT_EQ(13, getPidNucleusCharge(Al27));
  EXPECT_EQ(14, getPidNucleusCharge(Si28));
  EXPECT_EQ(14, getPidNucleusCharge(Si29));
  EXPECT_EQ(14, getPidNucleusCharge(Si30));
  EXPECT_EQ(15, getPidNucleusCharge(P31));
  EXPECT_EQ(16, getPidNucleusCharge(S32));
  EXPECT_EQ(16, getPidNucleusCharge(S33));
  EXPECT_EQ(16, getPidNucleusCharge(S34));
  EXPECT_EQ(17, getPidNucleusCharge(Cl35));
  EXPECT_EQ(18, getPidNucleusCharge(Ar36));
  EXPECT_EQ(17, getPidNucleusCharge(Cl37));
  EXPECT_EQ(18, getPidNucleusCharge(Ar38));
  EXPECT_EQ(19, getPidNucleusCharge(K39));
  EXPECT_EQ(20, getPidNucleusCharge(Ca40));
  EXPECT_EQ(19, getPidNucleusCharge(K41));
  EXPECT_EQ(20, getPidNucleusCharge(Ca42));
  EXPECT_EQ(20, getPidNucleusCharge(Ca43));
  EXPECT_EQ(20, getPidNucleusCharge(Ca44));
  EXPECT_EQ(21, getPidNucleusCharge(Sc45));
  EXPECT_EQ(22, getPidNucleusCharge(Ti46));
  EXPECT_EQ(22, getPidNucleusCharge(Ti47));
  EXPECT_EQ(22, getPidNucleusCharge(Ti48));
  EXPECT_EQ(22, getPidNucleusCharge(Ti49));
  EXPECT_EQ(24, getPidNucleusCharge(Cr50));
  EXPECT_EQ(23, getPidNucleusCharge(V51));
  EXPECT_EQ(24, getPidNucleusCharge(Cr52));
  EXPECT_EQ(24, getPidNucleusCharge(Cr53));
  EXPECT_EQ(26, getPidNucleusCharge(Fe54));
  EXPECT_EQ(25, getPidNucleusCharge(Mn55));
  EXPECT_EQ(26, getPidNucleusCharge(Fe56));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop