#ifndef SIMPROP_PHOTONFIELDS_PHOTONFIELD_H
#define SIMPROP_PHOTONFIELDS_PHOTONFIELD_H

namespace simprop {
namespace photonfield {

class AbstractField {
 public:
  AbstractField() {}
  virtual ~AbstractField() = default;

  virtual double getPhotonDensity(double ePhoton, double z = 0.) const { return 0.; }  // comoving
  std::string getFieldName() const { return m_fieldName; }
  void setEnhancementFactor(double f) { m_enhancementFactor = f; }

 protected:
  std::string m_fieldName = "AbstractPhotonField";
  double m_enhancementFactor = 1;
};

}  // namespace photonfield
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_PHOTONFIELD_H

// /**
//  @class TabularPhotonField
//  @brief Photon field decorator for tabulated photon fields.
//  This class reads photon field data from files;
//  The first file must be a list of photon energies [J], named fieldName_photonEnergy.txt
//  The second file must be a list of comoving photon field densities [1/m^3], named
//  fieldName_photonDensity.txt Optionally, a third file contains redshifts, named
//  fieldName_redshift.txt
//  */
// class TabularPhotonField : public PhotonField {
//  public:
//   TabularPhotonField(const std::string fieldName, const bool isRedshiftDependent = true);
//   double getPhotonDensity(double ePhoton, double z = 0.) const;
//   double getRedshiftScaling(double z) const;

//  protected:
//   void readPhotonEnergy(std::string filePath);
//   void readPhotonDensity(std::string filePath);
//   void readRedshift(std::string filePath);
//   void initRedshiftScaling();
//   void checkInputData() const;

//   std::vector<double> photonEnergies;
//   std::vector<double> photonDensity;
//   std::vector<double> redshifts;
//   std::vector<double> redshiftScalings;
// };

// TabularPhotonField::TabularPhotonField(std::string fieldName, bool isRedshiftDependent) {
// 	this->fieldName = fieldName;
// 	this->isRedshiftDependent = isRedshiftDependent;

// 	readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
// 	readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
// 	if (this->isRedshiftDependent)
// 		readRedshift(getDataPath("") + "Scaling/" + this->fieldName + "_redshift.txt");

// 	checkInputData();

// 	if (this->isRedshiftDependent)
// 		initRedshiftScaling();
// }

// double TabularPhotonField::getPhotonDensity(double ePhoton, double z) const {
// 	if (this->isRedshiftDependent) {
// 		return interpolate2d(ePhoton, z, this->photonEnergies, this->redshifts,
// this->photonDensity); 	} else { 		return interpolate(ePhoton,
// this->photonEnergies, this->photonDensity);
// 	}
// }

// double TabularPhotonField::getRedshiftScaling(double z) const {
// 	if (this->isRedshiftDependent) {
// 		if (z > this->redshifts.back()) {
// 			return 0.;
// 		} else if (z < this->redshifts.front()) {
// 			return 1.;
// 		} else {
// 			return interpolate(z, this->redshifts, this->redshiftScalings);
// 		}
// 	} else {
// 		return 1.;
// 	}
// }

// void TabularPhotonField::readPhotonEnergy(std::string filePath) {
// 	std::ifstream infile(filePath.c_str());
// 	if (!infile.good())
// 		throw std::runtime_error("TabularPhotonField::readPhotonEnergy: could not open " +
// filePath);

// 	std::string line;
// 	while (std::getline(infile, line)) {
// 		if (line.size() > 0)
// 			this->photonEnergies.push_back(std::stod(line));
// 	}
// 	infile.close();
// }

// void TabularPhotonField::readPhotonDensity(std::string filePath) {
// 	std::ifstream infile(filePath.c_str());
// 	if (!infile.good())
// 		throw std::runtime_error("TabularPhotonField::readPhotonDensity: could not open " +
// filePath);

// 	std::string line;
// 	while (std::getline(infile, line)) {
// 		if (line.size() > 0)
// 			this->photonDensity.push_back(std::stod(line));
// 	}
// 	infile.close();
// }

// void TabularPhotonField::readRedshift(std::string filePath) {
// 	std::ifstream infile(filePath.c_str());
// 	if (!infile.good())
// 		throw std::runtime_error("TabularPhotonField::initRedshift: could not open " +
// filePath);

// 	std::string line;
// 	while (std::getline(infile, line)) {
// 		if (line.size() > 0)
// 			this->redshifts.push_back(std::stod(line));
// 	}
// 	infile.close();
// }

// void TabularPhotonField::initRedshiftScaling() {
// 	double n0 = 0.;
// 	for (int i = 0; i < this->redshifts.size(); ++i) {
// 		double z = this->redshifts[i];
// 		double n = 0.;
// 		for (int j = 0; j < this->photonEnergies.size(); ++j) {
// 			double e = this->photonEnergies[j];
// 			if (z == 0.)
// 				n0 += getPhotonDensity(e, z);
// 			n += getPhotonDensity(e, z);
// 		}
// 		this->redshiftScalings.push_back(n / n0);
// 	}
// }

// void TabularPhotonField::checkInputData() const {
// 	if (this->isRedshiftDependent) {
// 		if (this->photonDensity.size() != this->photonEnergies.size() * this->
// redshifts.size()) 			throw
// std::runtime_error("TabularPhotonField::checkInputData: length of photon density input is unequal
// to length of photon energy input times length of redshift input"); 	} else { 		if
// (this->photonEnergies.size() != this->photonDensity.size()) 			throw
// std::runtime_error("TabularPhotonField::checkInputData: length of photon energy input is unequal
// to length of photon density input");
// 	}

// 	for (int i = 0; i < this->photonEnergies.size(); ++i) {
// 		double ePrevious = 0.;
// 		double e = this->photonEnergies[i];
// 		if (e <= 0.)
// 			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the
// photon energy input is not positive"); 		if (e <= ePrevious) throw
// std::runtime_error("TabularPhotonField::checkInputData: photon energy values are not strictly
// increasing"); 		ePrevious = e;
// 	}

// 	for (int i = 0; i < this->photonDensity.size(); ++i) {
// 		if (this->photonDensity[i] < 0.)
// 			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the
// photon density input is negative");
// 	}

// 	if (this->isRedshiftDependent) {
// 		if (this->redshifts[0] != 0.)
// 			throw std::runtime_error("TabularPhotonField::checkInputData: redshift input
// must start with zero");

// 		for (int i = 0; i < this->redshifts.size(); ++i) {
// 			double zPrevious = -1.;
// 			double z = this->redshifts[i];
// 			if (z < 0.)
// 				throw std::runtime_error("TabularPhotonField::checkInputData: a
// value in the redshift input is negative"); 			if (z <= zPrevious) throw
// std::runtime_error("TabularPhotonField::checkInputData: redshift values are not strictly
// increasing"); 			zPrevious = z;
// 		}

// 		for (int i = 0; i < this->redshiftScalings.size(); ++i) {
// 			double scalingFactor = this->redshiftScalings[i];
// 			if (scalingFactor <= 0.)
// 				throw std::runtime_error("TabularPhotonField::checkInputData:
// initRedshiftScaling has created a non-positive scaling factor");
// 		}
// 	}
// }
