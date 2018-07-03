#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "gp_matrix.hh"

void syntax_error() {
	fputs("Usage: ./output_text <in_file> <out_file>\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check for the correct number of command line arguments
	if(argc!=3) syntax_error();

	// Open input file, output contours, and delete temporary memory
	gp_matrix gp(argv[1]);

	FILE *fp=fopen(argv[2],"w");
	if(fp==NULL) {
		fputs("Error opening output file\n",stderr);
		return 1;
	}

	int i,j;
	for(j=0;j<gp.n;j++) {
		for(i=0;i<gp.m-1;i++) {
			fprintf(fp,"%g ",gp.f[i+gp.m*j]);
		}
		fprintf(fp,"%g\n",gp.f[gp.m-1+gp.m*j]);
	}
	fclose(fp);
}
