#include "sim_flatfold.hh"

#include <cmath>

int main() {
	sim_flatfold ff;

	// Create a diagonal normal vector
	double nx=2,ny=1,nr=sqrt(nx*nx+ny*ny);
	nx/=nr;ny/=nr;

	// Apply several folds
	ff.flatfold(nx,ny,0.1);
	ff.flatfold(0,1,0.1);
	ff.flatfold(0,1,-0.3);

	// Output positive and negative creases, plus the shape of the folded
	// sheet
	ff.crease_map("pos.dat",true);
	ff.crease_map("neg.dat",false);
	ff.output("ff.dat");
}
