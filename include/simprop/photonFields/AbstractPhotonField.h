#ifndef SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_PHOTONFIELD_H

namespace simprop {
namespace photonfield {

class AbstractField {
 public:
  AbstractField() {}
  virtual ~AbstractField() = default;

  virtual double getPhotonDensity(double ePhoton, double z = 0.) const = 0;  // comoving
  std::string getFieldName() const { return m_fieldName; }
  void setEnhancementFactor(double f) { m_enhancementFactor = f; }

 protected:
  std::string m_fieldName = "AbstractPhotonField";
  double m_enhancementFactor = 1;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
