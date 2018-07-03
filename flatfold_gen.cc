#include "sim_flatfold.hh"

#include <cstdlib>
#include <cmath>

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc<2||argc>5) {
		fputs("Syntax: ./flatfold_gen <num_folds> [<seed>] [<percent_radial>] [<sign>]\n\n"
		      "sign=0 for positive folds, sign=1 for folds with random sign\n",stderr);
		return 1;
	}

	// Check for a sensible number of folds
	int folds=atoi(argv[1]),i=0;
	if(folds<0||folds>100) {
		fputs("Fold number out of bounds\n",stderr);
		return 1;
	}

	// Check for optional arguments
	unsigned long seed=1;
	double frac=0.;
	bool rand_sign=false;
	if(argc>2) {

		// Read the random seed
		seed=atol(argv[2]);
		if(argc>3) {

			// Read the fraction of radial folds
			frac=0.01*atof(argv[3]);
			if(frac<0) frac=0;
			else if(frac>1) frac=1;

			// Read whether to use random signs on the folds
			if(argc>4) rand_sign=atoi(argv[4])==1;
		}
	}

	// Apply random folds, periodically reevaluating the bounding circle
	// for better efficiency in sampling folds
	sim_flatfold ff;
	ff.seed(seed);
	while(i<folds) {
		ff.random_fold(frac,rand_sign);
		if(++i%3==0) ff.compute_bounds();
	}

	// Output positive and negative creases, plus the shape of the folded
	// sheet
	ff.crease_map("pos.dat",true);
	ff.crease_map("neg.dat",false);
	ff.output("ff.dat");
}
