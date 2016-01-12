#ifndef TRIGGER_H
#define TRIGGER_H

#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>

namespace MOD {

   class Trigger {

      public:
         Trigger(std::string name, std::pair<int, int> prescales, bool fired);
         Trigger(std::istringstream & input_stream);
         Trigger();

         std::string name() const;
         std::string short_name() const;
         std::pair<int, int> prescale_pair() const;
         int prescale() const;
         bool fired() const;
         bool is_valid() const;
         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const Trigger&);

      private:
         std::string _name;
         std::pair<int, int> _prescales;
         bool _fired;
   };
}



#endif /* TRIGGER_H */