#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>

namespace MOD {

   class Property {

      public:
         Property(std::string name, int value);
         Property(std::string name, double value);
         Property(std::string name, std::string value);
         
         void value(std::string &) const;
         void value(int &) const;
         void value(double &) const;

         std::string name() const;

         std::string value_data_type() const;

         friend std::ostream& operator<< (std::ostream&, const Property&);

      private:
         std::string _name;
         std::string _data_type;

         int _int_value;
         double _double_value;
         std::string _string_value;

         std::string _value_data_type;
   };
}