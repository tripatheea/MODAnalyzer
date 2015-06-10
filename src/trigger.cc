#include "trigger.h"

using namespace std;

MOD::Trigger::Trigger(string name, pair<int, int> prescales, bool fired) : _name(name), _prescales(prescales), _fired(fired) {}

MOD::Trigger::Trigger(istringstream & input_stream) {
   string tag, name;
   bool fired;
   int prescale_1, prescale_2;

   input_stream >> tag >> name >> prescale_1 >> prescale_2 >> fired;

   _name = name;
   _prescales = make_pair( prescale_1, prescale_2 );
   _fired = fired;
}

MOD::Trigger::Trigger() : _fired(false) {}

string MOD::Trigger::name() const {
   return MOD::Trigger::_name;
}

pair<int, int> MOD::Trigger::prescale_pair() const {
   return MOD::Trigger::_prescales;
}

int MOD::Trigger::prescale() const {
   return _prescales.first * _prescales.second;
}

bool MOD::Trigger::fired() const {
   return MOD::Trigger::_fired;
}

bool MOD::Trigger::is_valid() const {
   return ( ! _name.empty());
}

string MOD::Trigger::make_string() const {
   stringstream ss;
   ss << "  Trig" 
        << setw(29) <<  _name 
        << setw(12) << _prescales.first 
        << setw(12) << _prescales.second 
        << setw(8) << _fired
        << endl;

   return ss.str();
}

string MOD::Trigger::make_header_string() const {
   stringstream ss;
   ss << "# Trig                         Name  Prescale_1  Prescale_2  Fired?" << endl;
   return ss.str();
}

namespace MOD {
  ostream& operator<< (ostream& os, const Trigger& trigger) {
     os << trigger.make_string();
     return os;
  }
}