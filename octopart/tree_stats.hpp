
#pragma once



struct tree_stats {
	int max_level;
	int nnodes;
	int nleaves;
	int nparts;
	void print() const;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & max_level;
		arc & nnodes;
		arc & nleaves;
		arc & nparts;
	}
};
