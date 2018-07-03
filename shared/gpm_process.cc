#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

void syntax_error() {
	fputs("Error processing command-line arguments\n",stderr);
	exit(2);
}

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc<3) {
		fputs("Usage: ./gpm_process [opts] <input_field> <output_field>\n\n"
		      "Carries out several operations on Gnuplot matrix fields\n\n"
		      "Options:\n"
		      "-c <phi_field>      Clean field using a level set\n"
		      "-p                  Project field to have zero mean\n"
		      "-r                  Reverse direction of cleaning\n"
		      "-n                  Use NaN to clean instead of zero\n"
		      "-k x0 x1 y0 y1      Crop to dimensions [x0:x1] [y0:y1]\n"
		      "-s                  Stepped output\n",stderr);
		return 1;
	}
	bool reverse=false,stepped=false,crop=false,project=false;
	char **ap,**ape=argv+argc-2,*clean=NULL;
	float val=0,x0,x1,y0,y1;

	// Process command-line options
	for(ap=argv+1;ap<ape;ap++) {
		if(strcmp("-c",*ap)==0) {
			clean=*(++ap);if(ap==ape) syntax_error();
		} else if(strcmp("-n",*ap)==0) val=std::numeric_limits<float>::quiet_NaN();
		else if(strcmp("-p",*ap)==0) project=true;
		else if(strcmp("-r",*ap)==0) reverse=true;
		else if(strcmp("-k",*ap)==0) {
			if(ap+4>=ape) syntax_error(); // Possibly slighty dodgy
			x0=atof(*(++ap));x1=atof(*(++ap));
			y0=atof(*(++ap));y1=atof(*(++ap));
			crop=true;
		} else if(strcmp("-s",*ap)==0) stepped=true;
		else {
			printf("gpm_process: unrecognized option \"%s\"\n",*ap);
			return 1;
		}
	}

	// Open input file and levelset file
	gp_matrix fld(*(ap++));

	// Project field to have zero mean if requested
	if(project) fld.zero_project();

	// Clean the field and output
	if(clean!=NULL) {
		int i;
		gp_matrix phi(clean);
		if(fld.m==phi.m&&fld.n==phi.n) {
			for(i=0;i<fld.mn;i++) if(reverse?phi.f[i]<0:phi.f[i]>0) fld.f[i]=val;
		} else if(fld.m==phi.m+1&&fld.n==phi.n+1) {
			int j,ij;
			float ph;
			for(ij=j=0;j<fld.n;j++) for(i=0;i<fld.m;i++,ij++) {
				ph=phi.field_lp_int(i-0.5,j-0.5);
				if(reverse?ph<0:ph>0) fld.f[ij]=val;
			}
		} else {
			fprintf(stderr,"Field dimensions (%d,%d) don't agree with phi dimensions (%d,%d)\n"
					,fld.m,fld.n,phi.m,phi.n);
			return 1;
		}
	}

	// Crop the field and output the stepped field if necessary
	if(crop) fld.crop(x0,x1,y0,y1);
	stepped?fld.output_stepped(*ap):fld.output(*ap);
}
