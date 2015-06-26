#include "condition.h"

using namespace std;

MOD::Condition::Condition(int lumi_block, double avg_inst_lumi, int npv) : _lumi_block(lumi_block), _avg_inst_lumi(avg_inst_lumi), _npv(npv) {}

MOD::Condition::Condition(istringstream & input_stream) {
   string tag;
   int lumi_block, npv;
   double avg_inst_lumi;

   input_stream >> tag >> lumi_block >> avg_inst_lumi >> npv;

   _lumi_block = lumi_block;
   _avg_inst_lumi = avg_inst_lumi;
   _npv = npv;
}

MOD::Condition::Condition() {}

int MOD::Condition::lumi_block() const {
  return _lumi_block;
}

double MOD::Condition::avg_inst_lumi() const {
  return _avg_inst_lumi;
}

int MOD::Condition::npv() const {
  return _npv;
}


string MOD::Condition::make_string() const {
   stringstream ss;

    ss << "  Cond "
       << setw(9) << _lumi_block
       << setw(12) << _avg_inst_lumi
       << setw(4) << _npv
       << endl;   

   return ss.str();
}

string MOD::Condition::make_header_string() const {
   stringstream ss;
   ss << "# Cond LumiBlock AvgInstLumi NPV" << endl;
   return ss.str();
}

namespace MOD {
  ostream& operator<< (ostream& os, const Condition& cond) {
     os << cond.make_string();
     return os;
  }
}