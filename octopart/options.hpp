#pragma once
#include <string>

class options {
public:
	std::string config_file;
	bool dust_only;
	bool first_order_space;
	bool first_order_time;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & config_file;
		arc & dust_only;
		arc & first_order_space;
		arc & first_order_time;
	}

	bool process_options(int argc, char *argv[]);
};
