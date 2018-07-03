#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "gp_matrix.hh"

void syntax_error() {
	fputs("Usage: ./make_contour [-p <phi>] <in_file> <out_file> {<height> <height> ...}\n"
	      "       ./make_contour [-p <phi>] <in_file> <out_file> r <start> <incr> <n>\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check for the correct number of command line arguments
	if(argc<3) syntax_error();

	// Check if a level set file is given
	int q;
	if(strcmp(argv[1],"-p")==0) {q=2;if(argc<5) syntax_error();}
	else q=0;

	// Process the command-line arguments to obtain the heights
	float *hm;int hc;
	if(argc==q+3) {
		hm=new float[1];*hm=0;hc=1;
	} else {
		if(strcmp(argv[q+3],"r")==0) {
			if(argc!=7+q) syntax_error();
			hc=atoi(argv[6+q]);
			if(hc<=0||hc>16777216) {
				fprintf(stderr,"Value of n=%d invalid\n",hc);
				return 1;
			}
			hm=new float[hc];
			double a=atof(argv[4+q]),b=atof(argv[5+q]);
			for(int i=0;i<hc;i++) hm[i]=a+i*b;
		} else {
			hm=new float[hc=argc-3-q];
			for(int i=3+q;i<argc;i++) hm[i-3-q]=atof(argv[i]);
		}
	}

	// Open input file, output contours, and delete temporary memory
	gp_matrix gp(argv[q+1]);
	if(q==2) {
		gp_matrix phi(argv[2]);
		gp.output_contours(argv[4],hm,hc,phi);
	} else gp.output_contours(argv[2],hm,hc);
	delete [] hm;
}
