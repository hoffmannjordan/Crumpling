include ../config.mk

# Lists of files to be built
objs=mat3.o mat.o gp_matrix.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mat_test gpm_process gpm_metrics make_contour bitmap_field colorbar gpm_extrema output_text

all: lib $(execs)

include Makefile.dep

depend:
	$(cxx) -MM $(src) >Makefile.dep

clean:
	rm $(objs) $(execs)

gp_matrix.o: gp_matrix.cc
	$(im_cxx) $(im_cflags) -c $<

%.o: %.cc
	$(cxx) $(cflags) -c $<

lib: libmatvec.a libgpmtx.a

libmatvec.a: mat.o mat3.o
	rm -f libmatvec.a
	ar rs libmatvec.a $^

libgpmtx.a: gp_matrix.o
	rm -f libgpmtx.a
	ar rs libgpmtx.a $^

output_text: output_text.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

gpm_process: gpm_process.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

gpm_metrics: gpm_metrics.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

gpm_extrema: gpm_extrema.cc libgpmtx.a
	$(im_cxx) $(cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

make_contour: make_contour.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

bitmap_field: bitmap_field.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

colorbar: colorbar.cc libgpmtx.a
	$(im_cxx) $(im_cflags) -o $@ $< -L. -lgpmtx $(im_lflags)

mat_test: mat_test.cc libmatvec.a
	$(cxx) $(cflags) -o $@ $< -L. -lmatvec

.PHONY: clean depend lib
