

class MODNTilde {

	public:
		MODNTilde(double cone_radius, int pt_cut);
		double calculate_n_tilde(MODEvent * event);

	private:
		double _cone_radius;
		int _pt_cut;
};