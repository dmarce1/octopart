#pragma once

#include <octopart/dim.hpp>

class nindex {
	int i;
public:
	operator int();
	nindex begin();
	nindex end();
	int get_dim(int d) const;
	void set_dim(int d, int v);
};
