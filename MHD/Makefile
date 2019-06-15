ifndef config
config = magdif.in
endif
magdifplot := ../RUN/PLOTTING/magdifplot.py
mesh := ../PRELOAD/inputformaxwell.msh

Bn_diffs := $(sort $(wildcard Bn_diff_*.dat))
Bn_sums := $(patsubst Bn_diff_%.dat, Bn_%.dat, $(Bn_diffs))

extract_proj = sed -e 's/ \+/\t/g; s/^\t//g' $(1) | cut -f 9-10 > $(2)
poloidal_modes = cd ../RUN && ./result_spectrum.x $(1) > $(2) && cd ../MHD

.PHONY: run plots
.ONESHELL: conv_sum.dat conv_diff.dat

run: $(config) $(mesh)
	./magdif_test.sh $(config)

plots: Bmn_psi.dat Bmn_vac_psi.dat conv_sum.dat conv_diff.dat \
	$(magdifplot) $(config) $(mesh)
	python $(magdifplot) . $(config) $(mesh)

Bmn_psi.dat: plot_Bn.dat
	$(call extract_proj, $<, ../RUN/Bn_psi.dat)
	$(call poloidal_modes, Bn_psi.dat, $@)
	mv ../RUN/$@ ./

Bmn_vac_psi.dat: plot_Bn_vac.dat
	$(call extract_proj, $<, ../RUN/Bn_vac_psi.dat)
	$(call poloidal_modes, Bn_vac_psi.dat, $@)
	mv ../RUN/$@ ./

conv_sum.dat: $(Bn_sums)
	rm -f $@
	cd ../FEM
	FreeFem++-nw L2int.edp ../MHD/Bn_vac.dat ../MHD/$@
	for file in $^ ; do FreeFem++-nw L2int.edp ../MHD/$$file ../MHD/$@ ; done
	cd ../MHD

conv_diff.dat: $(Bn_diffs)
	rm -f $@
	cd ../FEM
	for file in $^ ; do FreeFem++-nw L2int.edp ../MHD/$$file ../MHD/$@ ; done
	cd ../MHD
