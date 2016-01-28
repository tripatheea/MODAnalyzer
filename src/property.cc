#include "property.h"

using namespace std;

MOD::Property::Property(std::string name, int value) : _name(name), _int_value(value), _value_data_type("int") {

}

MOD::Property::Property(std::string name, double value) : _name(name), _double_value(value), _value_data_type("double") {

}

MOD::Property::Property(std::string name, std::string value) : _name(name), _string_value(value), _value_data_type("string") {

}

MOD::Property::Property(std::string name, long value) : _name(name), _long_value(value), _value_data_type("long") {

}

void MOD::Property::value(int & value) const {
  value = _int_value;
}

void MOD::Property::value(double & value) const {
  value = _double_value;
}

void MOD::Property::value(long & value) const {
  value = _long_value;
}

void MOD::Property::value(string & value) const {
  value = _string_value;
}

std::string MOD::Property::name() const {
  return _name;
}

string MOD::Property::value_data_type() const {
  return _value_data_type;
}

namespace MOD {
  
  ostream& operator<< (ostream& os, const Property& property) {
    
    if (property.value_data_type() == "string") {
      std::string v;
      property.value(v);
      os << v;
    }
    else if (property.value_data_type() == "int") {
      int v;
      property.value(v);
      os << v;
    }
    else if (property.value_data_type() == "double") {
      double v;
      property.value(v);
      os << v;
    }
    else if (property.value_data_type() == "long") {
      long v;
      property.value(v);
      os << v;
    }
    
    return os;
  }

}