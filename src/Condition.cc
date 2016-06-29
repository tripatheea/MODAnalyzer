#include "Condition.h"

using namespace std;

MOD::Condition::Condition(istringstream & input_stream) {
   string tag;
   int run_number, event_number, lumi_block, npv;
   double average_instantaneous_lumi, integrated_delivered_lumi, integrated_recorded_lumi;
   long timestamp, ms_offset;
   bool valid_lumi;

   input_stream >> tag >> run_number >> event_number >> lumi_block >> valid_lumi >> integrated_delivered_lumi >> integrated_recorded_lumi >> average_instantaneous_lumi >> npv >> timestamp >> ms_offset;

   _run_number = run_number;
   _event_number = event_number;
   _lumi_block = lumi_block;
   _valid_lumi = valid_lumi;
   _integrated_delivered_lumi = integrated_delivered_lumi;
   _integrated_recorded_lumi = integrated_recorded_lumi;
   _average_instantaneous_lumi = average_instantaneous_lumi;
   _npv = npv;
   _timestamp = timestamp;
   _ms_offset = ms_offset;
}

MOD::Condition::Condition() {}

const int MOD::Condition::lumi_block() const {
  return _lumi_block;
}

const double MOD::Condition::average_instantaneous_lumi() const {
  return _average_instantaneous_lumi;
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

const bool MOD::Condition::valid_lumi() const {
  return _valid_lumi;
}

const double MOD::Condition::integrated_delivered_lumi() const {
  return _integrated_delivered_lumi;
}

const double MOD::Condition::integrated_recorded_lumi() const {
  return _integrated_recorded_lumi;
}

const long MOD::Condition::time() const {
  long recorded_time = _timestamp; 
  // recorded_time = recorded_time << 32;
  // recorded_time += _ms_offset;

  return recorded_time;
}

string MOD::Condition::make_string() const {
   stringstream ss;
    ss << "    Cond"
       << setw(16) << _run_number
       << setw(16) << _event_number
       << setw(16) << _lumi_block
       << setw(16) << _valid_lumi
       << setw(16) << _integrated_delivered_lumi
       << setw(16) << _integrated_recorded_lumi
       << setw(16) << _average_instantaneous_lumi
       << setw(16) << _npv
       << setw(16) << _timestamp
       << setw(16) << _ms_offset
       << endl;   

   return ss.str();
}

string MOD::Condition::make_header_string() const {
   stringstream ss;
   ss << "#   Cond          RunNum        EventNum       LumiBlock       validLumi     intgDelLumi     intgRecLumi     AvgInstLumi             NPV       timestamp        msOffset" << endl;
   return ss.str();
}

namespace MOD {
  ostream& operator<< (ostream& os, const Condition& cond) {
     os << cond.make_string();
     return os;
  }
}