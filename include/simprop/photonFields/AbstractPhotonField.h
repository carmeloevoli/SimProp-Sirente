#ifndef SIMPROP_PHOTONFIELDS_ABSTRACTPHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_ABSTRACTPHOTONFIELD_H

namespace simprop {
namespace photonfield {

class AbstractPhotonField {
 public:
  AbstractPhotonField(std::string fieldName, double ePhotonMin, double ePhotonMax)
      : m_fieldName(fieldName), m_ePhotonMin(ePhotonMin), m_ePhotonMax(ePhotonMax) {}
  virtual ~AbstractPhotonField() = default;

  virtual double getPhotonDensity(double ePhoton, double z = 0.) const = 0;  // comoving
  virtual double I_gamma(double ePhoton, double z = 0.) const = 0;

  std::string getFieldName() const { return m_fieldName; }
  double getMinPhotonEnergy() const { return m_ePhotonMin; }
  double getMaxPhotonEnergy() const { return m_ePhotonMax; }

 protected:
  std::string m_fieldName = "AbstractPhotonField";
  double m_ePhotonMin;
  double m_ePhotonMax;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
