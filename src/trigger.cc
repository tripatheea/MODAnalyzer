#include <iostream>


using namespace std;

class Trigger {

	public:
		Trigger(string name, pair<int, int> prescales, bool fired);
		Trigger();

		string name();
		pair<int, int> prescales();
		bool fired();
		bool is_valid();

	private:
		string name_;
		bool fired_;
		pair<int, int> prescales_;		
};

Trigger::Trigger(string name, pair<int, int> prescales, bool fired) : name_(name), prescales_(prescales), fired_(fired) {}
Trigger::Trigger() {}

string Trigger::name() {
	return Trigger::name_;
}

pair<int, int> Trigger::prescales() {
	return Trigger::prescales_;
}

bool Trigger::fired() {
	return Trigger::fired_;
}

bool Trigger::is_valid() {
	return ( ! name_.empty());
}