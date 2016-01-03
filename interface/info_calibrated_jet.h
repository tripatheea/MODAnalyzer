#ifndef InfoCalibratedJet_H
#define InfoCalibratedJet_H

#include <iostream>

#include "fastjet/ClusterSequence.hh"

namespace MOD {

	class InfoCalibratedJet : public fastjet::PseudoJet::UserInfoBase {

		public:
			enum JetQualityLevels_t { UNDETERMINED = -1, FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3 };
			
			InfoCalibratedJet(std::string tag, double JEC, double area, int number_of_constituents, int charged_multiplicity, double neutral_hadron_fraction, double neutral_em_fraction, double charged_hadron_fraction, double charged_em_fraction, double eta);
			InfoCalibratedJet(std::string tag);
			InfoCalibratedJet(std::string tag, double JEC);

			const std::string tag() const;
			const double JEC() const;
						
			const double area() const;

			const int number_of_constituents() const;
			const int charged_multiplicity() const;
			const double neutral_hadron_fraction() const;
			const double neutral_em_fraction() const;
			const double charged_hadron_fraction() const;
			const double charged_em_fraction() const;

			const JetQualityLevels_t jet_quality() const;

			void set_jet_quality_level();

			const std::string header() const;
			
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

			
			JetQualityLevels_t _jet_quality;

	};
}


#endif /* InfoCalibratedJet_H */