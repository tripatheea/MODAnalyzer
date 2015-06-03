#include "../interface/trigger.h"

using namespace std;

MODTrigger::MODTrigger(string name, pair<int, int> prescales, bool fired) : _name(name), _prescales(prescales), _fired(fired) {}

MODTrigger::MODTrigger(string input_string) {
	vector<string> components = split(input_string);

	_name = components[1];
	_prescales = make_pair( stoi(components[2]), stoi(components[3]) );
	_fired = (stoi(components[4]) == 1);
}

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

vector<string> MODTrigger::split(string const &input) { 
	istringstream buffer(input);
	vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
	return ret;
}