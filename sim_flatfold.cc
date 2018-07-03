#include "sim_flatfold.hh"

#include <cmath>

/** Initializes the flatfold class to be a square sheet covering [-1,1]^2. */
sim_flatfold::sim_flatfold() : rng(gsl_rng_alloc(gsl_rng_taus)) {
	f.push_back(new facet(-1,1,-1,1));
	cx=2;cy=0;crsq=2.;cr=sqrt(2.);
}

/** Initializes the flatfold class to be a rectangular sheet.
 * \param[in] (ax,bx) the x-dimensions of the sheet.
 * \param[in] (ay,by) the y-dimensions of the sheet. */
sim_flatfold::sim_flatfold(double ax,double bx,double ay,double by) :
	rng(gsl_rng_alloc(gsl_rng_taus)) {
	f.push_back(new facet(ax,bx,ay,by));
	cx=0.5*(bx+ax);cy=0.5*(by+ay);
	double lx=bx-ax,ly=by-ay;
	crsq=0.25*(lx*ly+ly*ly);
	cr=sqrt(crsq);
}

/** The class destructor frees the dynamically allocated facets. */
sim_flatfold::~sim_flatfold() {
	for(unsigned int i=0;i<f.size();i++) delete f[i];
	gsl_rng_free(rng);
}

/** Computes a bounding circle to the current sheet, by first finding the
 * center of mass, and then finding the maximum distance to a vertex. */
void sim_flatfold::compute_bounds() {
	cx=0;cy=0;
	double sa=0,rsq;

	// Find the center of mass
	for(unsigned int i=0;i<f.size();i++)
		f[i]->weight_contrib(cx,cy,sa);
	cx/=sa;cy/=sa;

	// Compute the maximum distance of a vertex to the center of mass. For
	// quick lookups, both the distance and the square distance are stored.
	crsq=0.;
	for(unsigned int i=0;i<f.size();i++)
		if((rsq=f[i]->max_rad_sq(cx,cy))>crsq) crsq=rsq;
	cr=sqrt(crsq);
}

/** Applies a random flat fold to the sheet, choosing a random angle and
 * displacement.
 * \param[in] rand_sign whether to choose a random sign for the fold or not.
 * \return Whether the fold intersected the sheet or not. */
bool sim_flatfold::random_flatfold(bool rand_sign) {
	double th=2*M_PI*gsl_rng_uniform(rng);
	double nx=sin(th),ny=cos(th);
	return flatfold(nx,ny,cx*nx+cy*ny+cr*(-1+2*gsl_rng_uniform(rng)),
			rand_sign?random_sign():1);
}

/** Applies a random radial fold to the sheet.
 * \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_radial_fold(bool rand_sign) {
	int k=0;
	double x,y;
	while(k<sim_flatfold_max_attempts) {

		// Find a point inside the bounding circle
		x=-1+2*gsl_rng_uniform(rng);
		y=-1+2*gsl_rng_uniform(rng);
		if(x*x+y*y>crsq) continue;

		// Scale the point and check if it is inside the folded sheet.
		// If so, apply a random radial fold.
		x=cx+cr*x;
		y=cy+cr*y;
		if(point_inside(x,y)) {
			radial_fold(x,y,2*M_PI*gsl_rng_uniform(rng),
				    M_PI*gsl_rng_uniform_pos(rng),
				    M_PI*gsl_rng_uniform_pos(rng),
				    rand_sign?random_sign():1);
			return;
		}
		k++;
	}
	fputs("Too many attempts to find a point inside the sheet\n",stderr);
	exit(1);
}

/** Applies a flat fold to the sheet.
 * \param[in] (nx,ny) the normal vector of the fold to apply.
 * \param[in] di the displacement of the fold to apply.
 * \param[in] fsign the sign of the fold.
 * \return Whether the line given actually folded the sheet or not. */
bool sim_flatfold::flatfold(double nx,double ny,double di,int fsign) {
	unsigned int i=0,n=f.size();

	// Check to see if the fold intersects any facet. If not, then skip it
	for(i=0;i<n;i++) if(f[i]->intersects(nx,ny,di)) break;
	if(i==n) return false;

	// If the code reaches here, then we know that the fold intersects a
	// facet. Perform the folding by splitting the facet into two, and
	// mirroring one about the line of the fold.
	for(i=0;i<n;i++) {
		facet *nf=new facet(*f[i]);
		bool p1=f[i]->c.nplane(nx,ny,di,fsign);
		if(nf->c.nplane(-nx,-ny,-di,0)) {
			nf->mirror(nx,ny,di);

			// Update the bounding circle, since the mirrored facet
			// might extend further than before
			update_bounding_circle(nf);
			if(p1) f.push_back(nf);
			else {
				delete f[i];
				f[i]=nf;
			}
		} else {
			delete nf;
			if(!p1) {
				fputs("Error: both facets removed by fold\n",stderr);
				exit(1);
			}
		}
	}
	return true;
}

/** Applies a radial fold to the sheet.
 * \param[in] (x,y) the focus point of the radial fold.
 * \param[in] ro the overall rotation of the fold position.
 * \param[in] (al,be) the first two angles making up the radial fold.
 * \param[in] fsign the sign of the three folds that share the same sign. */
void sim_flatfold::radial_fold(double x,double y,double ro,double al,double be,int fsign) {
	unsigned int i=0,j,n=f.size();
	facet *qf[4];

	// Compute the normals to the four planes making up this radial fold
	double nxa=cos(ro),nya=sin(ro),nxb=cos(ro+al),nyb=sin(ro+al),
	       nxc=cos(ro+al+be),nyc=sin(ro+al+be),nxd=cos(ro+be+M_PI),nyd=sin(ro+be+M_PI),
	       da=nxa*x+nya*y,db=nxb*x+nyb*y,dc=nxc*x+nyc*y,dd=nxd*x+nyd*y;

	for(;i<n;i++) {
		j=0;
		facet *nbe=new facet(*f[i]),*nga=new facet(*f[i]),*nde=new facet(*f[i]);

		// Compute the first facet
		if(f[i]->c.nplane(-nxa,-nya,-da,fsign)&&f[i]->c.nplane(nxb,nyb,db,fsign)) qf[j++]=f[i];
		else delete f[i];

		// Compute the second facet
		if(nbe->c.nplane(-nxb,-nyb,-db,0)&&nbe->c.nplane(nxc,nyc,dc,0)) {
			nbe->mirror(nxb,nyb,db);
			update_bounding_circle(nbe);
			qf[j++]=nbe;
		} else delete nbe;

		// Compute the third facet
		if(nga->c.nplane(-nxc,-nyc,-dc,fsign)&&nga->c.nplane(nxd,nyd,dd,0)) {
			nga->mirror(nxc,nyc,dc);
			nga->mirror(nxb,nyb,db);
			update_bounding_circle(nga);
			qf[j++]=nga;
		} else delete nga;

		// Compute the fourth facet
		if(nde->c.nplane(-nxd,-nyd,-dd,-fsign)&&nde->c.nplane(nxa,nya,da,0)) {
			nde->mirror(nxa,nya,da);
			update_bounding_circle(nde);
			qf[j++]=nde;
		} else delete nde;

		// At least once facet should be present. If not, print an
		// error message
		if(j==0) {
			fputs("Error: all four facets removed by radial fold\n",stderr);
			exit(1);
		}

		// Store the available facets. One can be placed at the
		// location of the original facet. Others need to be added.
		f[i]=*qf;
		while(j>1) f.push_back(qf[--j]);
	}
}

/** Applies a random fold to the sheet.
 * \param[in] frac the probability of the fold being a radial fold.
 * \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold(double frac,bool rand_sign) {
	if(gsl_rng_uniform(rng)>frac) {
		for(int k=0;k<sim_flatfold_max_attempts;k++)
			if(random_flatfold(rand_sign)) return;
		fputs("Too many flatfold attempts in random_fold\n",stderr);
		exit(1);
	} else random_radial_fold(rand_sign);
}

/** Updates the bounding circle to ensure in encompasses a given facet.
 * \param[in] nf a pointer to the facet. */
void sim_flatfold::update_bounding_circle(facet *nf) {
	double nrsq=nf->max_rad_sq(cx,cy);
	if(nrsq>crsq) {crsq=nrsq;cr=sqrt(crsq);}
}

/** Checks to see if a given position is inside any facet.
 * \param[in] (x,y) the position to consider.
 * \return Whether the position is inside or not. */
bool sim_flatfold::point_inside(double x,double y) {
	for(unsigned int i=0;i<f.size();i++)
		if(f[i]->point_inside(x,y)) return true;
	return false;
}

/** Outputs the positions of creases in the original coordinate system.
 * \param[in] fp the file handle to write to.
 * \param[in] positive whether to write positive or negative creases. */
void sim_flatfold::crease_map(FILE *fp,bool positive) {
	for(unsigned int i=0;i<f.size();i++)
		f[i]->unfolded_crease(fp,positive);
}

/** Outputs the total crease mileage.
 * \param[out] pos the positive crease mileage.
 * \param[out] neg the negative crease mileage. */
void sim_flatfold::crease_mileage(double &pos,double &neg) {
	pos=neg=0;
	for(unsigned int i=0;i<f.size();i++)
		f[i]->crease_mileage(pos,neg);
}

/** Outputs the facet positions in the current coordinate system.
 * \param[in] fp the file handle to write to. */
void sim_flatfold::output(FILE *fp) {
	for(unsigned int i=0;i<f.size();i++) f[i]->output(fp);
}
