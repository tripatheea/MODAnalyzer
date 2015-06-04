#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>

#include "fastjet/ClusterSequence.hh"


class MODParticle {

   public:
      MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, std::string trigger_type);
      MODParticle(std::istringstream & input_stream);
      MODParticle();

      fastjet::PseudoJet pseudojet() const;
      int pdgId() const;
      double mass() const;
      std::string make_string() const;
      std::string make_header_string() const;

      friend std::ostream& operator<< (std::ostream&, const MODParticle&);

   private:
      std::string _trigger_type;
      fastjet::PseudoJet _pseudojet;
      double _mass;
      int _pdgId; 
};