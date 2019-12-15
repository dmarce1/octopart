#pragma once
#include <octopart/real.hpp>
#include <string>

class options {
public:
	std::string config_file;
	bool dust_only;
	bool first_order_space;
	bool first_order_time;
	bool fpe;
	bool periodic;
	bool reflecting;
	int parts_per_node;
	int problem_size;
	double fgamma;
	double tmax;
	std::string problem;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & config_file;
		arc & dust_only;
		arc & first_order_space;
		arc & first_order_time;
		arc & fpe;
		arc & parts_per_node;
		arc & periodic;
		arc & reflecting;
		arc & problem_size;
		arc & fgamma;
		arc & tmax;
		arc & problem;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
