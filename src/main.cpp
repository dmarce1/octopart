/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include <hpx/hpx_init.hpp>

#include <octopart/tree.hpp>

int hpx_main(int argc, char *argv[]) {
	std::vector<particle>parts;
	auto root = hpx::new_<tree>(hpx::find_here(), std::move(parts), null_range()).get();
	typename tree::set_self_action f;
	f(root, root);
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	hpx::init(argc, argv);
}