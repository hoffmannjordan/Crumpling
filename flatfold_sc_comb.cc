#include <cstdio>
#include <cmath>

const int max_fold=66;

int main(int argc,char **argv) {

	// Print syntax message if no files are provided
	if(argc<2) {
		fputs("Syntax: ./flatfold_sc_comb <dat1> <dat2> ...\n",stderr);
		return 1;
	}

	// Initialize accumulators
	double sl[6*max_fold],*sll=sl+max_fold,
	       *sm=sll+max_fold,*smm=sm+max_fold,
	       *sd=smm+max_fold,*sdd=sd+max_fold,
	       sf[6*max_fold];
	int ttrials=0,trials,k;
	for(k=0;k<6*max_fold;k++) sl[k]=0;

	// Loop over all available data files
	for(int i=1;i<argc;i++) {

		// Open and read file contents, checking for errors or missing data
		FILE *fp=fopen(argv[i],"rb");
		if(fp==NULL) {
			fputs("Can't open output file\n",stderr);
			return 1;
		}
		if(fread(&trials,sizeof(int),1,fp)!=1) {
			fputs("File input error\n",stderr);
			return 1;
		}
		printf("# File %s: %d trials\n",argv[i],trials);
		if(fread(sf,sizeof(double),6*max_fold,fp)!=(size_t) 6*max_fold) {
			fputs("File input error\n",stderr);
			return 1;
		}

		// Add the contributions from this file to the accumulators
		ttrials+=trials;
		for(k=0;k<6*max_fold;k++) sl[k]+=sf[k];
		fclose(fp);
	}

	// Print the combined statistics
	double ifac=1./ttrials;
	for(k=0;k<max_fold;k++) {
		double ml=sl[k]*ifac,mm=sm[k]*ifac,md=sd[k]*ifac;
		printf("%d %g %g %g %g %g %g\n",k,ml,sqrt(sll[k]*ifac-ml*ml),
		       mm,sqrt(smm[k]*ifac-mm*mm),md,sqrt(sdd[k]*ifac-md*md));
	}
}
