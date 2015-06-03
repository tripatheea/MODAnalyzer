#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>


class MODTrigger {

	public:
		MODTrigger(std::string name, std::pair<int, int> prescales, bool fired);
		MODTrigger(std::string input_string);
		MODTrigger();

		std::string name() const;
		std::pair<int, int> prescale_pair() const;
		int prescale() const;
		bool fired() const;
		bool is_valid() const;
		std::string make_string() const;
		std::string make_header_string() const;

		friend std::ostream& operator<< (std::ostream&, const MODTrigger&);

	private:
		std::string _name;
		bool _fired;
		std::pair<int, int> _prescales;

		std::vector<std::string> split(std::string const &input);
};