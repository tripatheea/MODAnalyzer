#include <iostream>


using namespace std;

class MODTrigger {

	public:
		MODTrigger(string name, pair<int, int> prescales, bool fired);
		MODTrigger();

		string name();
		pair<int, int> prescale_pair();
		int prescale();
		bool fired();
		bool is_valid();

	private:
		string _name;
		bool _fired = false;
		pair<int, int> _prescales;		
};

MODTrigger::MODTrigger(string name, pair<int, int> prescales, bool fired) : _name(name), _prescales(prescales), _fired(fired) {}
MODTrigger::MODTrigger() {}

string MODTrigger::name() {
	return MODTrigger::_name;
}

pair<int, int> MODTrigger::prescale_pair() {
	return MODTrigger::_prescales;
}

int MODTrigger::prescale() {
	return _prescales.first * _prescales.second;
}

bool MODTrigger::fired() {
	return MODTrigger::_fired;
}

bool MODTrigger::is_valid() {
	return ( ! _name.empty());
}