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
	if(argc<2) {
		fputs("Usage: ./gpm_extrema [opts] <input_field>\n\n"
		      "Calculates the extremal value of a Gnuplot matrix fields\n\n"
		      "Options:\n"
		      "-c <phi_field>      Clean field using a level set\n"
		      "-r                  Reverse direction of cleaning\n",stderr);
		return 1;
	}
	bool reverse=false;
	char **ap,**ape=argv+argc-1,*clean=NULL;
	float vmin,vmax;

	// Process command-line options
	for(ap=argv+1;ap<ape;ap++) {
		if(strcmp("-c",*ap)==0) {
			clean=*(++ap);if(ap==ape) syntax_error();
		} else if(strcmp("-r",*ap)==0) reverse=true;
		else {
			printf("gpm_extrema: unrecognized option \"%s\"\n",*ap);
			return 1;
		}
	}

	// Open input file and levelset file
	gp_matrix fld(*(ap++));

	// Clean the field and output
	if(clean!=NULL) {
		int i;bool start=true;
		gp_matrix phi(clean);
		if(fld.m==phi.m&&fld.n==phi.n) {
			for(i=0;i<fld.mn;i++) if(reverse?phi.f[i]>0:phi.f[i]<0) {
				if(start) {
					vmin=vmax=fld.f[i];start=false;
				} else {
					if(fld.f[i]<vmin) vmin=fld.f[i];
					else if (fld.f[i]>vmax) vmax=fld.f[i];
				}
			}
		} else if(fld.m==phi.m+1&&fld.n==phi.n+1) {
			int j,ij;
			float ph;
			for(ij=j=0;j<fld.n;j++) for(i=0;i<fld.m;i++,ij++) {
				ph=phi.field_lp_int(i-0.5,j-0.5);
				if(reverse?ph>0:ph<0) {
					if(start) {
						vmin=vmax=fld.f[ij];start=false;
					} else {
						if(fld.f[ij]<vmin) vmin=fld.f[ij];
						else if (fld.f[ij]>vmax) vmax=fld.f[ij];
					}
				}
			}
		} else {
			fprintf(stderr,"Field dimensions (%d,%d) don't agree with phi dimensions (%d,%d)\n"
					,fld.m,fld.n,phi.m,phi.n);
			return 1;
		}
		if(start) {
			fputs("No valid gridpoints found\n",stderr);
			return 1;
		}
	} else {
		vmin=vmax=fld.f[0];
		for(float *fp=fld.f+1,*fe=fp+fld.mn;fp<fe;fp++) {
			if(*fp<vmin) {vmin=*fp;}
			else if(*fp>vmax) vmax=*fp;
		}
	}

	// Output the extremal values
	printf("%g %g\n",vmin,vmax);
}
