#ifndef InfoPFC_H
#define InfoPFC_H

#include <iostream>

#include "fastjet/ClusterSequence.hh"

namespace MOD {

	class InfoPFC : public fastjet::PseudoJet::UserInfoBase {

		public:
			InfoPFC(int pdgId, std::string tag);
			
			const int pdgId() const;
			const std::string tag() const;
			
		private:
			int _pdgId;
			std::string _tag;

	};
}


#endif /* InfoPFC_H */