#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>

namespace MOD {

   class Condition {

      public:
         Condition(int lumi_block, double avg_inst_lumi, int npv);
         Condition(std::istringstream & input_stream);
         Condition();

         int lumi_block() const;
         double avg_inst_lumi() const;
         int npv() const;

         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const Condition&);

      private:
         int _lumi_block;
         double _avg_inst_lumi;
         int _npv;
   };
}