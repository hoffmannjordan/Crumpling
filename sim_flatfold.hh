#ifndef SIM_FLATFOLD_HH
#define SIM_FLATFOLD_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gsl/gsl_rng.h"

#include "facet.hh"

/** The maximum number of attempts for cases where a random trial result
 * satisfying conditions must be found. */
const int sim_flatfold_max_attempts=65536;

class sim_flatfold {
	public:
		/** A vector of polygonal facets the make up the sheet. */
		std::vector<facet*> f;
		sim_flatfold();
		sim_flatfold(double ax,double bx,double ay,double by);
		~sim_flatfold();
		bool random_flatfold(bool rand_sign=false);
		void random_radial_fold(bool rand_sign=false);
		void random_fold(double frac,bool rand_sign=false);
		bool flatfold(double nx,double ny,double di,int fsign=1);
		void radial_fold(double x,double y,double ro,double al,double be,int fsign=1);
		void crease_map(FILE *fp,bool positive);
		inline void crease_map(const char* filename,bool positive) {
			FILE *fp=fopen(filename,"w");
			if(fp==NULL) {
				fputs("Can't open output file\n",stderr);
				exit(1);
			}
			crease_map(fp,positive);
			fclose(fp);
		}
		void crease_mileage(double &pos,double &neg);
		void output(FILE *fp);
		inline void output(const char* filename) {
			FILE *fp=fopen(filename,"w");
			if(fp==NULL) {
				fputs("Can't open output file\n",stderr);
				exit(1);
			}
			output(fp);
			fclose(fp);
		}
		void compute_bounds();
		/** Seeds the random number generator.
		 * \param[in] seed the seed. */
		inline void seed(unsigned long seed) {
			gsl_rng_set(rng,seed);
		}
		bool point_inside(double x,double y);
	private:
		/** Picks a random sign for a fold.
		 * \return The sign. */
		inline int random_sign() {
			return gsl_rng_get(rng)&1?-1:1;
		}
		void update_bounding_circle(facet *f);
		/** The x-coordinate of the approximate center of the folded
		 * sheet. */
		double cx;
		/** The y-coordinate of the approximate center of the folded
		 * sheet. */
		double cy;
		/** The maximum square radius of a vertex from the approximate
		 * center. */
		double crsq;
		/** The maximum radius of a vertex from the approximate center.
		 */
		double cr;
		/** A GSL random number generator. */
		gsl_rng* const rng;
};

#endif
