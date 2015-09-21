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
         Condition(int run_number, int event_number, int lumi_block, double avg_inst_lumi, int npv);
         Condition(std::istringstream & input_stream);
         Condition();

         const int lumi_block() const;
         const double avg_inst_lumi() const;
         const int npv() const;

         const int event_number() const;
         const int run_number() const;

         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const Condition&);

      private:
         int _run_number;
         int _event_number;

         int _lumi_block;
         double _avg_inst_lumi;
         int _npv;
         
   };
}


#endif /* CONDITION_H */