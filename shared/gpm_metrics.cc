#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc!=2) {
		fputs("Usage: ./gpm_metrics <input_field>\n\n"
		      "Reads in a gnuplot matrix file and calculates the grid spacings\n",stderr);
		return 1;
	}

	// Read in the Gnuplot file and compute the grid spacings, assuming the
	// grid is equally spaced. Use the maximum extent of the grid to
	// minimize roundoff errors.
	gp_matrix fld(argv[1]);
	int &m=fld.m,&n=fld.n;
	printf("%g %g\n",(fld.x[m-1]-*fld.x)/(m-1),(fld.y[n-1]-*fld.y)/(n-1));
}
