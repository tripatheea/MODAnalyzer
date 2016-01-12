#include "condition.h"

using namespace std;

MOD::Condition::Condition(int run_number, int event_number, int lumi_block, double avg_inst_lumi, int npv) : _run_number(run_number), _event_number(event_number), _lumi_block(lumi_block), _avg_inst_lumi(avg_inst_lumi), _npv(npv) {}

MOD::Condition::Condition(istringstream & input_stream) {
   string tag;
   int run_number, event_number, lumi_block, npv;
   double avg_inst_lumi;

   input_stream >> tag >> run_number >> event_number >> lumi_block >> avg_inst_lumi >> npv;

   _run_number = run_number;
   _event_number = event_number;
   _lumi_block = lumi_block;
   _avg_inst_lumi = avg_inst_lumi;
   _npv = npv;
}

MOD::Condition::Condition() {}

const int MOD::Condition::lumi_block() const {
  return _lumi_block;
}

const double MOD::Condition::avg_inst_lumi() const {
  return _avg_inst_lumi;
}

const int MOD::Condition::npv() const {
  return _npv;
}

const int MOD::Condition::event_number() const {
  return _event_number;
}

const int MOD::Condition::run_number() const {
  return _run_number;
}

string MOD::Condition::make_string() const {
   stringstream ss;
    ss << "    Cond"
       << setw(16) << _run_number
       << setw(16) << _event_number
       << setw(16) << _lumi_block
       << setw(16) << _avg_inst_lumi
       << setw(16) << _npv
       << endl;   

   return ss.str();
}

string MOD::Condition::make_header_string() const {
   stringstream ss;
   ss << "#   Cond          RunNum        EventNum       LumiBlock     AvgInstLumi             NPV" << endl;
   return ss.str();
}

namespace MOD {
  ostream& operator<< (ostream& os, const Condition& cond) {
     os << cond.make_string();
     return os;
  }
}