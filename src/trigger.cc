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
		string make_string();
		string header();

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

string MODTrigger::make_string() {
	stringstream ss;
	ss << "trig" 
		  << setw(16) << _name 
		  << setw(15) << _prescales.first 
		  << setw(20) << _prescales.second 
		  << setw(17) << _fired
		  << endl;
		  
	return ss.str();
}

string MODTrigger::header() {
	stringstream ss;
	ss << "#Trig          Name          Prescale_1          Prescale_2          Fired?" << endl;
	return ss.str();
}