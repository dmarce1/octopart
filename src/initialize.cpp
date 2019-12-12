#include <octopart/initialize.hpp>


void drift(particle& p) {
	p.m = 1.0;
	p.e = 0.0;
	const auto r = abs(p.x);
	p.u = -p.x;
}


init_func_type get_initialization_function(const std::string& name) {
	if( name == "drift") {
		return drift;
	} else {
		printf( "Error: Initialization function %s is not known\n", name.c_str());
		abort();
	}
}
