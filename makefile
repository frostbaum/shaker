FC=gfortran
FCFLAGS+= -O2
p_NAME := shakerv2.1

.PHONY: all, clean, distclean

all: $(p_NAME)

$(p_NAME): v3d_func_rep.o rng.o structure_c.o geo_data_c.o shaker.o main.o
	$(FC) $(FCFLAGS) v3d_func_rep.o rng.o structure_c.o geo_data_c.o shaker.o main.o -o $(p_NAME)

v3d_func_rep.o v3d_func_rep.mod: 3d-vectorOPs/v3d_func_rep.f95
	$(FC) $(FCFLAGS) -c 3d-vectorOPs/v3d_func_rep.f95
	@ touch v3d_func_rep.mod

rng.o rng.mod: rng/rng.f95 v3d_func_rep.mod
	$(FC) $(FCFLAGS) -c rng/rng.f95
	@ touch rng_c.mod

structure_c.o structure_c.mod: structure/structure_c.f95 v3d_func_rep.mod
	$(FC) $(FCFLAGS) -c structure/structure_c.f95
	@ touch structure_c.mod

geo_data_c.o geo_data_c.mod: geo_data_c.f95
	$(FC) $(FCFLAGS) -c geo_data_c.f95
	@ touch geo_data_c.mod

shaker.o shaker.mod: shaker.f95 structure_c.mod geo_data_c.mod rng.mod
	$(FC) $(FCFLAGS) -c shaker.f95
	@ touch shaker.mod

main.o: main.f95 shaker.mod
	$(FC) $(FCFLAGS) -c main.f95

clean:
	@- $(RM) $(p_NAME)
	@- $(RM) *.o *.mod

distclean: clean
