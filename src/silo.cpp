/*
 * silo.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#include <cstdio>
#include <silo.h>
#include <octopart/options.hpp>
#include <octopart/delaunay.hpp>
#include <octopart/math.hpp>
#include <octopart/particle.hpp>
#include <octopart/math.hpp>

int main(int argc, char *argv[]) {
	static auto opts = options::get();
	FILE *fp = fopen(argv[1], "rb");
	if (fp == NULL) {
		printf("File not found\n");
	} else {
		std::vector<particle> parts;
		particle p;
		fixed_real tm;
		std::string filename = argv[1];
		filename += ".silo";
		int cnt = 0;
		cnt += fread(&opts.fgamma, sizeof(real), 1, fp);
		cnt += fread(&tm, sizeof(real), 1, fp);
		opts.set(opts);
		while (p.read(fp)) {
			parts.push_back(p);
		}
		fclose(fp);
		auto optlist = DBMakeOptlist(2);
		float ftime = (float) (double) tm;
		double dtime = tm;
		DBAddOption(optlist, DBOPT_TIME, &ftime);
		DBAddOption(optlist, DBOPT_DTIME, &dtime);
		DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Meshless",
		DB_HDF5);
#if( NDIM == 1)
		const int shapetypes[1] = {DB_ZONETYPE_BEAM};
	#elif( NDIM ==2)
		const int shapetypes[1] = { DB_ZONETYPE_TRIANGLE };
#else
		const int shapetypes[1] = { DB_ZONETYPE_TET };
#endif
		const int shapesize[1] = { NDIM + 1 };
		auto delaunay_regions = compute_delaunay_regions(parts);
		const int nzones = delaunay_regions.size();
		const int nnodes = parts.size();
		const int shapecount[1] = { nzones };
		std::vector<int> node_list;
		node_list.reserve(NDIM * delaunay_regions.size());
		for (auto &r : delaunay_regions) {
#if(NDIM==3)
			const auto a = parts[r[1]].x - parts[r[0]].x;
			const auto b = parts[r[2]].x - parts[r[0]].x;
			const auto c = parts[r[3]].x - parts[r[0]].x;
			if (triple_product(c, a, b) > 0.0) {
				std::swap(r[0], r[1]);
			}
#endif
			for (int n = 0; n < NDIM + 1; n++) {
				node_list.push_back(r[n]);
			}
		}
		real *coords[NDIM];
		char *coordnames[NDIM];
		for (int dim = 0; dim < NDIM; dim++) {
			coords[dim] = new real[nnodes];
			coordnames[dim] = new char[2];
			coordnames[dim][0] = 'x' + dim;
			coordnames[dim][1] = '\0';
			for (int i = 0; i < nnodes; i++) {
				coords[dim][i] = parts[i].x[dim];
			}
		}
		DBPutZonelist2(db, "zonelist", nzones, NDIM, node_list.data(), node_list.size(), 0, 0, 0, shapetypes, shapesize, shapecount, 1,
		optlist);
		DBPutUcdmesh(db, "mesh", NDIM, coordnames, coords, nnodes, nzones, "zonelist", NULL, DB_DOUBLE, optlist);
		std::vector<real> Nc;
		std::vector<real> h;
		std::vector<real> t;
		std::vector<real> dt;
		std::vector<real> rho;
		std::vector<real> pre;
		std::array<std::vector<real>, NDIM> vel;
		std::array<std::vector<real>, NDIM> g;
		std::array<std::vector<real>, NDIM> dpre_dx;
		t.reserve(nnodes);
		dt.reserve(nnodes);
		rho.reserve(nnodes);
		pre.reserve(nnodes);
			for (int dim = 0; dim < NDIM; dim++) {
			vel[dim].reserve(nnodes);
			dpre_dx[dim].reserve(nnodes);
			g[dim].reserve(nnodes);
		}
		for (const auto &pi : parts) {
			const auto Q = pi.Q;
			t.push_back(double(pi.t));
			dt.push_back(double(pi.dt));
			rho.push_back(pi.W.rho);
			pre.push_back(pi.W.p);
			h.push_back(pi.h);
			std::array<vect, NDIM> E;
			Nc.push_back(pi.Nc);
			for (int dim = 0; dim < NDIM; dim++) {
				vel[dim].push_back(pi.W.v[dim]);
				dpre_dx[dim].push_back(pi.dW[dim].p);
				g[dim].push_back(pi.g[dim]);
			}
		}
		DBPutUcdvar1(db, "t", "mesh", t.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		DBPutUcdvar1(db, "dt", "mesh", dt.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		DBPutUcdvar1(db, "h", "mesh", h.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		DBPutUcdvar1(db, "Nc", "mesh", Nc.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		DBPutUcdvar1(db, "rho", "mesh", rho.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		DBPutUcdvar1(db, "p", "mesh", pre.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "dpre_d" + char('x' + char(dim));
			DBPutUcdvar1(db, nm.c_str(), "mesh", dpre_dx[dim].data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "v_" + char('x' + char(dim));
			DBPutUcdvar1(db, nm.c_str(), "mesh", vel[dim].data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "g_" + char('x' + char(dim));
			DBPutUcdvar1(db, nm.c_str(), "mesh", g[dim].data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			delete[] coordnames[dim];
			delete[] coords[dim];
		}
		DBFreeOptlist(optlist);
		DBClose(db);

	}

}
