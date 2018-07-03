#include "facet.hh"

#include <cstdlib>

using namespace voro;

/** Initializes a polygonal facet to be a complete rectangular sheet.
 * \param[in] (ax,bx) the x dimensions of the sheet.
 * \param[in] (ay,by) the y dimensions of the sheet. */
facet::facet(double ax,double bx,double ay,double by) :
	flipped(false), mtrans(1), vtrans(0) {
	c.init_base(0.5*ax,0.5*bx,0.5*ay,0.5*by);
	for(int i=0;i<4;i++) c.ne[i]=0;
}

/** Initializes a polygonal facet to be a copy of another one.
 * \param[in] f the facet to be copied. */
facet::facet(facet &f) : flipped(f.flipped),
	mtrans(f.mtrans.a,f.mtrans.b,f.mtrans.c,f.mtrans.d),
	vtrans(f.vtrans.x,f.vtrans.y), c(f.c) {}

/** Computes the contribution of this facet to the center of mass
 * of the folded sheet.
 * \param[in,out] (sx,sy) the accumulators for the first moments,
 *			  which have the contribution from this
 *			  facet added to them.
 * \param[in,out] sa the accumulator for the area, which will have the
 *	             contribution from this facet added to it. */
void facet::weight_contrib(double &sx,double &sy,double &sa) {
	double cx,cy,ca;

	// Compute centroid and area of this facet. Apply a scaling of two due
	// to the rescalings done by the Voronoi library. Take the absolute
	// value of the area since some facets may be flipped in orientation.
	c.centroid(cx,cy);cx*=2;cy*=2;
	ca=fabs(c.area())*4;
	sx+=cx*ca;sy+=cy*ca;
	sa+=ca;
}

/** Computes the maximum radius squared of a vertex to a given location.
 * \param[in] (cx,cy) the location to use.
 * \return The maxium radius squared. */
double facet::max_rad_sq(double cx,double cy) {
	double dx,dy,mrsq=0,rsq,*pp=c.pts;
	for(int i=0;i<c.p;i++) {
		dx=cx-*(pp++);
		dy=cy-*(pp++);
		rsq=dx*dx+dy*dy;
		if(rsq>mrsq) mrsq=rsq;
	}
	return mrsq;
}

/** Determines if this facet intersects a given plane.
 * \param[in] (nx,ny) the normal vector of the plane.
 * \param[in] di the displacement of the plane. */
bool facet::intersects(double nx,double ny,double di) {
	return c.plane_intersects(nx,ny,di)==c.plane_intersects(-nx,-ny,-di);
}

/** Mirrors the facet about a given plane.
 * \param[in] (nx,ny) the normal vector of the plane.
 * \param[in] di the displacement of the plane. */
void facet::mirror(double nx,double ny,double di) {
	flipped=!flipped;
	sym_mat m(1-2*nx*nx,-2*nx*ny,1-2*ny*ny);
	vec v(2*di*nx,2*di*ny);
	for(int i=0;i<c.p;i++) {

		// Mirror the vertex position
		vec w(c.pts[2*i],c.pts[2*i+1]);
		w=m*w+v;
		c.pts[2*i]=w.x;
		c.pts[2*i+1]=w.y;

		// Flip the sign of any crease information
		c.ne[i]=-c.ne[i];
	}

	// Update transformation matrix for this facet
	mtrans=m*mtrans;
	vtrans=m*vtrans+v;
}

/** Outputs the positions of creases from this facet in the original
 * coordinate system.
 * \param[in] fp the file handle to write to.
 * \param[in] positive whether to write positive or negative creases. */
void facet::unfolded_crease(FILE *fp,bool positive) {
	if(c.p==0) {
		fputs("Error: zero vertex facet\n",stderr);
		exit(1);
	}

	// Loop around the vertices, and check for specific value 1 or -1 in
	// the neighbor information
	int nev=positive!=flipped?1:-1,k=0;
	vec v1,v2;
	do {
		if(c.ne[k]==nev) {
			v1=transform(k);
			v2=transform(k=c.ed[2*k]);
			fprintf(fp,"%g %g\n%g %g\n\n",v1.x,v1.y,v2.x,v2.y);
		} else k=c.ed[2*k];
	} while(k!=0);
}

/** Checks whether a given position is inside the facet.
 * \param[in] (x,y) the position.
 * \return Whether the position is inside or not. */
bool facet::point_inside(double x,double y) {
	int k=0,l,w=0;
	double a,b;
	bool left=*c.pts<x,nleft;

	// Compute the winding number by counting the number of times the facet
	// edges cross a line from (x,y) pointing upward
	do {
		l=c.ed[2*k];
		nleft=c.pts[2*l]<x;
		if(nleft!=left) {
			a=c.pts[2*l+1]*(x-c.pts[2*k])+c.pts[2*k+1]*(c.pts[2*l]-x);
			b=y*(c.pts[2*l]-c.pts[2*k]);
			if(left) {if(a<b) w++;}
			else {if(a>b) w--;}
		}
		left=nleft;
		k=l;
	} while(k!=0);

	// Flipped facets will have the winding number reversed. Correct for
	// this here.
	if(flipped) w=-w;

	// The winding number should be either 0 (if the point is outside) or 1
	// (if the point is inside). If it's not one of these cases, then print
	// an error message.
	if(w<0) {
		fprintf(stderr,"Negative winding number of %d for (%g,%g)\n",w,x,y);
		exit(1);
	} else if(w>1) {
		fprintf(stderr,"Winding number of %d for (%g,%g)\n",w,x,y);
		exit(1);
	}
	return w>0;
}

/** Adds contributions to the crease mileage from this facet.
 * \param[in,out] pos the cumulative positive crease mileage.
 * \param[in,out] neg the cumulative negative crease mileage. */
void facet::crease_mileage(double &pos,double &neg) {
	if(c.p==0) {
		fputs("Error: zero vertex facet\n",stderr);
		exit(1);
	}

	// Loop around the vertices, and check for specific value 1 or -1 in
	// the neighbor information
	int k=0;
	vec v1,v2;
	do {
		if(c.ne[k]!=0) {
			v1=transform(k);
			bool p=(c.ne[k]==1)!=flipped;
			(p?pos:neg)+=sqrt(mod_sq(v1-transform(k=c.ed[2*k])));
		} else k=c.ed[2*k];
	} while(k!=0);
}

/** Outputs the facet position in the current coordinate system.
 * \param[in] fp the file handle to write to. */
void facet::output(FILE *fp) {
	if(c.p==0) return;
	int k=0;

	// Loop around the vertices and print them
	do {
		fprintf(fp,"%g %g\n",c.pts[2*k],c.pts[2*k+1]);
		k=c.ed[2*k];
	} while (k!=0);

	// Print the first one again, to make a complete loop
	fprintf(fp,"%g %g\n\n",c.pts[0],c.pts[1]);
}

/** Transforms a given vertex from the current coordinate system back to the
 * original one.
 * \param[in] k the index of the vertex to consider.
 * \return The transformed position vector. */
vec facet::transform(int k) {
	vec v(c.pts[2*k],c.pts[2*k+1]);
	v=mtrans.transpose()*(v-vtrans);
	return v;
}
