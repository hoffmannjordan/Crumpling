#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "bi_interp.hh"

// The radius to integrate to
const double rad=1.;

// The number of Simpson's rule evaluation points
const int sam=100;

// The fraction of function values to save to the output file
const int every=4;

// The integration step size
const double h=2*rad/sam;

int main() {
	double *u=new double[49],cth,sth,x,y,f,sum,sum2;

	// Create some test data to interpolate in a checkboard pattern
	for(int j=0;j<7;j++) for(int i=0;i<7;i++)
		u[i+7*j]=(i+j)&1;

	// Construct the bicubic interpolation class
	bicubic_interp bi(u,7,7,0,6,0,6);

	// Open files to output the diagnostic information
	FILE *fp=fopen("btest1.out","w"),
	     *fp2=fopen("btest2.out","w");
	if(fp==NULL||fp2==NULL) {
		fputs("Error opening output files\n",stderr);
		return 1;
	}

	// Loop over a range of angles. Output the field data on rays with each
	// angle. Integrate with Simpson's rule and compare to the exact answer
	// based on quadrature.
	for(double th=0;th<0.995*M_PI;th+=0.01*M_PI) {
		cth=cos(th);sth=sin(th);
		sum=0.;
		for(int i=0;i<=sam;i++) {

			// Evaluate the function, print it if needed, and
			// compute its contribution to the Simpson's rule
			// calculation
			x=3+(i*h-rad)*cth;y=3+(i*h-rad)*sth;
			f=bi.f(x,y);
			if(i%every==0) fprintf(fp,"%g %g %.10g\n",x,y,bi.f(x,y));
			sum+=i&1?2.*f:(i==0||i==sam?0.5*f:f);
		}

		// Finish the Simpson's rule calculation and compare to the
		// exact result
		sum*=2*h/3.;
		sum2=bi.line_integral(3-rad*cth,3-rad*sth,3+rad*cth,3+rad*sth);
		fprintf(fp2,"%.10g %.12g %.12g %.12g\n",th,sum,sum2,sum-sum2);
	}

	// Close the files and free the dynamically allocated memory
	delete [] u;
}
