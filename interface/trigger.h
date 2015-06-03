#include <iostream>
#include <vector>

class MODTrigger {

	public:
		MODTrigger(std::string name, std::pair<int, int> prescales, bool fired);
		MODTrigger(std::string input_string);
		MODTrigger();

		const std::string name() const;
		const std::pair<int, int> prescale_pair() const;
		const int prescale() const;
		const bool fired() const;
		const bool is_valid() const;
		const std::string make_string();
		const std::string header() const;

	private:
		std::string _name;
		bool _fired = false;
		std::pair<int, int> _prescales;

		std::vector<std::string> split(std::string const &input);
};