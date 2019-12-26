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
	bool gravity;
	bool periodic;
	bool reflecting;
	int parts_per_node;
	int problem_size;
	double fgamma;
	double kep_eps;
	double theta;
	double tmax;
	double cfl;
	double output_freq;
	std::string problem;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & cfl;
		arc & kep_eps;
		arc & checkpoint;
		arc & config_file;
		arc & dust_only;
		arc & first_order_space;
		arc & first_order_time;
		arc & fpe;
		arc & gravity;
		arc & parts_per_node;
		arc & periodic;
		arc & reflecting;
		arc & problem_size;
		arc & fgamma;
		arc & theta;
		arc & tmax;
		arc & problem;
		arc & output_freq;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
