#pragma once
#include <octopart/real.hpp>
#include <string>

class options {
public:
	std::string config_file;
	std::string checkpoint;
	bool dust_only;
	bool first_order_space;
	bool first_order_time;
	bool fpe;
	bool global_time;
	bool gravity;
	bool lloyd_correct;
	bool x_periodic;
	bool y_periodic;
	bool z_periodic;
	bool x_reflecting;
	bool y_reflecting;
	bool z_reflecting;
	int parts_per_node;
	int problem_size;
	double fgamma;
	double kep_eps;
	double theta;
	double tmax;
	double cfl;
	double output_freq;
	double grid_size;
	std::string problem;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & lloyd_correct;
		arc & cfl;
		arc & kep_eps;
		arc & checkpoint;
		arc & config_file;
		arc & dust_only;
		arc & first_order_space;
		arc & first_order_time;
		arc & fpe;
		arc & global_time;
		arc & gravity;
		arc & parts_per_node;
		arc & x_periodic;
		arc & y_periodic;
		arc & z_periodic;
		arc & x_reflecting;
		arc & y_reflecting;
		arc & z_reflecting;
		arc & problem_size;
		arc & fgamma;
		arc & theta;
		arc & tmax;
		arc & problem;
		arc & output_freq;
		arc & grid_size;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
