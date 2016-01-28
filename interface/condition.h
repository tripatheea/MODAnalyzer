#ifndef CONDITION_H
#define CONDITION_H

#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>

namespace MOD {

   class Condition {

      public:
         Condition(std::istringstream & input_stream);
         Condition();

         const int lumi_block() const;
         const double average_instantaneous_lumi() const;
         const int npv() const;

         const bool valid_lumi() const;
         const double integrated_delivered_lumi() const;
         const double integrated_recorded_lumi() const;
         const long time() const;


         const int event_number() const;
         const int run_number() const;

         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const Condition&);

      private:
         int _run_number;
         int _event_number;

         int _lumi_block;
         double _average_instantaneous_lumi;
         int _npv;

         bool _valid_lumi;
         double _integrated_delivered_lumi;
         double _integrated_recorded_lumi;
         long _timestamp;
         long _ms_offset;

         
   };
}


#endif /* CONDITION_H */