#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

#ifdef HAS_IMAGEMAGICK
#include "Magick++.h"
using namespace Magick;
#endif

int main(int argc,char **argv) {

#ifdef HAS_IMAGEMAGICK
	// Initialize ImageMagick - this is required on some installations
	// before using any ImageMagick functions
	InitializeMagick(*argv);
#endif

	// Check for the correct number of command-line arguments
	if(argc<4||argc>5) {
		fputs("Usage: ./colorbar <output_png> <width> <height> [<palette>]\n",stderr);
		return 1;
	}

	// Check the dimensions are sensible
	int wid=atoi(argv[2]),hei=atoi(argv[3]),ptype=0;
	if(wid<=0||wid>16777216||hei<=0||hei>16777216) {
		fputs("Image dimensions out of bounds\n",stderr);
		return 1;
	}

	// Check the palette number, if specified
	if(argc==5) {
		ptype=atoi(argv[4]);
		if(ptype<0||ptype>4) {
			fputs("Palette number out of range\n",stderr);
			return 1;
		}
	}

	// Create and output the color bar
	gp_matrix_bitmap_colorbar(argv[1],wid,hei,ptype);
}
