include ./config.mk

lflags=`gsl-config --libs`
iflags=`gsl-config --cflags` -I./shared -I./voro++_2d

objs=sim_flat_fold.o facet.o
#objs=parab_func.o bi_interp.o bi_image.o facet.o sim_flatfold.o
src=$(patsubst %.o,%.cc,$(objs))
#execs=gsl_min_example matio_example im_test bic_test ls_test \
      radon_test imaging_model flatfold_test flatfold_gen \
      flatfold_scale flatfold_sc_comb im_extract
execs=flatfold_test flatfold_gen
all:
	$(MAKE) -C ./shared
	$(MAKE) -C ./voro++_2d
	$(MAKE) executables

executables: $(execs)

include Makefile.dep

depend:
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

#im_test: im_test.cc bi_interp.o bi_image.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lmatio

#im_extract: im_extract.cc bi_interp.o bi_image.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lmatio

#ls_test: ls_test.cc bi_interp.o bi_image.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lmatio

#radon_test: radon_test.cc bi_interp.o bi_image.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lmatio

flatfold_test: flatfold_test.cc sim_flatfold.o facet.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -L./shared -lmatvec -L./voro++_2d -lvoro++_2d

flatfold_gen: flatfold_gen.cc sim_flatfold.o facet.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -L./shared -lmatvec -L./voro++_2d -lvoro++_2d

#flatfold_scale: flatfold_scale.cc sim_flatfold.o facet.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -L./shared -lmatvec -L./voro++_2d -lvoro++_2d

flatfold_sc_comb: flatfold_sc_comb.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^

#bic_test: bic_test.cc bi_interp.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

#gsl_min_example: gsl_min_example.cc parab_func.o
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

#imaging_model: imaging_model.cc
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

#matio_example: matio_example.cc
#	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lmatio

.PHONY: clean depend executables
