#ifndef InfoCalibratedJet_H
#define InfoCalibratedJet_H

#include <iostream>

#include "fastjet/ClusterSequence.hh"

namespace MOD {

	class InfoCalibratedJet : public fastjet::PseudoJet::UserInfoBase {

		public:
			InfoCalibratedJet(std::string tag, double JEC, double area, int number_of_constituents, int charged_multiplicity, double neutral_hadron_fraction, double neutral_em_fraction, double charged_hadron_fraction, double charged_em_fraction, double eta);
			
			const std::string tag() const;
			const double JEC() const;
			const int number_of_constituents() const;
			const int jet_quality();
			
		private:
			std::string _tag;

			double _JEC;
			double _area;

			int _number_of_constituents;
			int _charged_multiplicity;

			double _neutral_hadron_fraction;
			double _neutral_em_fraction;
			double _charged_hadron_fraction;
			double _charged_em_fraction;

			double _eta;

			enum JetQualityLevels_t { UNDETERMINED = -1, FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3 };
			JetQualityLevels_t _jet_quality;

	};
}


#endif /* InfoCalibratedJet_H */