#ifndef HELPERS_H
#define HELPERS_H


#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>
#include <stdexcept>
#include <cmath>


#include "fastjet/ClusterSequence.hh"
#include "InfoPFC.h"
#include "InfoCalibratedJet.h"






namespace MOD {
 
  std::vector<fastjet::PseudoJet> filter_charged(std::vector<fastjet::PseudoJet>);

}


#endif /* HELPERS_H */