#ifndef SIMPROP_PHOTONFIELDS_ABSTRACTPHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_ABSTRACTPHOTONFIELD_H

namespace simprop {
namespace photonfield {

class AbstractPhotonField {
 public:
  AbstractPhotonField() {}
  virtual ~AbstractPhotonField() = default;

  virtual double density(double ePhoton, double z = 0.) const = 0;
  virtual double I_gamma(double ePhoton, double z = 0.) const = 0;

  double getMinPhotonEnergy() const { return m_ePhotonMin; }
  double getMaxPhotonEnergy() const { return m_ePhotonMax; }

 protected:
  double m_ePhotonMin = 0.;
  double m_ePhotonMax = 0.;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
