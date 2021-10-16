#ifndef SIMPROP_PID_H
#define SIMPROP_PID_H

#include <cassert>

namespace simprop {

class PID {
 private:
  int m_id = -1;
  std::string m_name = "none";

 public:
  explicit PID() {}

  explicit PID(int Z, int A, std::string name) : m_name(name) {
    assert(Z <= A);
    m_id = idNucleus(Z, A);
  }

  explicit PID(int elementaryId, std::string name) : m_name(name) {
    assert(elementaryId < 10);
    m_id = elementaryId;
  }

  int idNucleus(int Z, int A) { return 1000000000 + 10 * Z + 10000 * A; }
  bool isGood() const { return (m_id > 0); }
  bool isNucleus() const { return (m_id > 1000000000); }  // TODO check for ap this
  bool operator==(const PID& other) const { return m_id == other.m_id; }
  bool operator!=(const PID& other) const { return m_id != other.m_id; }
  bool operator<(const PID& other) const { return m_id < other.m_id; }
  bool operator>(const PID& other) const { return m_id > other.m_id; }

  friend std::ostream& operator<<(std::ostream& stream, const PID& pid) {
    stream << pid.m_name;
    return stream;
  }
};

static const PID photon(0, "photon");
static const PID neutrino_e(1, "nu_e");
static const PID antineutrino_e(2, "nubar_e");
static const PID neutrino_mu(3, "nu_mu");
static const PID antineutrino_mu(4, "nubar_mu");
static const PID electron(5, "e-");
static const PID positron(6, "e+");

static const PID neutron(0, 1, "neutron");
static const PID proton(1, 1, "proton");
static const PID antiproton(-1, 1, "antiproton");
static const PID He4(2, 4, "He4");
static const PID C12(6, 12, "C12");
static const PID N14(7, 14, "N14");
static const PID O16(8, 16, "O16");
static const PID Fe56(26, 56, "Fe56");

}  // namespace simprop

#endif  // SIMPROP_PID_H