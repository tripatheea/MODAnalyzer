#include <iostream>
#include <vector>

class MODTrigger {

	public:
		MODTrigger(std::string name, std::pair<int, int> prescales, bool fired);
		MODTrigger(std::string input_string);
		MODTrigger();

		std::string name();
		std::pair<int, int> prescale_pair();
		int prescale();
		bool fired();
		bool is_valid();
		std::string make_string();
		std::string header();

	private:
		std::string _name;
		bool _fired = false;
		std::pair<int, int> _prescales;

		std::vector<std::string> split(std::string const &input);
};