#ifndef FACET_HH
#define FACET_HH

#include "vec.hh"
#include "mat.hh"

#include "cell_2d.hh"

#include <cstdio>

struct facet {
	/** Whether the facet is current flipped over from its starting
	 * configuration. */
	bool flipped;
	/** The (orthogonal) matrix transformation from the original facet
	 * position to the current position. */
	mat mtrans;
	/** The vector displacement from the original facet position to the
	 * current position. */
	vec vtrans;
	/** A Voronoi cell, which represents this facet as a convex polyhedron.
	 */
	voro::voronoicell_neighbor_2d c;
	facet(double ax,double bx,double ay,double by);
	facet(facet &f);
	void weight_contrib(double &sx,double &sy,double &sa);
	bool intersects(double nx,double ny,double di);
	void mirror(double nx,double ny,double di);
	bool point_inside(double x,double y);
	double max_rad_sq(double cx,double cy);
	void unfolded_crease(FILE *fp,bool positive);
	void crease_mileage(double &pos,double &neg);
	void output(FILE *fp);
	vec transform(int k);
};

#endif
