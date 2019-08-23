ifndef config
config = magdif.in
endif
magdifplot := ../RUN/PLOTTING/magdifplot.py
mesh := ../PRELOAD/inputformaxwell.msh

Bn_diffs := $(sort $(wildcard Bn_diff_*.dat))
Bn_sums := $(patsubst Bn_diff_%.dat, Bn_%.dat, $(Bn_diffs))

extract_proj = sed -e 's/ \+/\t/g; s/^\t//g' $(1) | cut -f 9-10 > $(2)
poloidal_modes = cd ../RUN && ./result_spectrum.x $(1) > $(2) && cd ../MHD

.PHONY: run plots clean
.ONESHELL: conv_sum.dat conv_diff.dat

run: $(config) $(mesh)
	./magdif_test.sh $(config)

plots: Bmn_psi.dat Bmn_vac_psi.dat conv_sum.dat conv_diff.dat \
	$(magdifplot) $(config) $(mesh)
	python $(magdifplot) . $(config) $(mesh)

clean:
	rm -f *.pdf Bn*.dat plot*.dat *.out presn*.dat currn*.dat fluxvar.dat fort.333

Bmn_psi.dat: plot_Bn.dat
	$(call extract_proj, $<, Bn_psi.dat)
	$(call poloidal_modes, ../MHD/Bn_psi.dat, ../MHD/$@)
	rm Bn_psi.dat

Bmn_vac_psi.dat: plot_Bn_vac.dat
	$(call extract_proj, $<, Bn_vac_psi.dat)
	$(call poloidal_modes, ../MHD/Bn_vac_psi.dat, ../MHD/$@)
	rm Bn_vac_psi.dat

conv_sum.dat: $(Bn_sums)
	rm -f $@
	cd ../FEM
	FreeFem++-nw L2int.edp ../MHD/Bn_vac.dat ../MHD/$@ 1>/tmp/freefem.out 2>&1
	for file in $^
	do
	FreeFem++-nw L2int.edp ../MHD/$$file ../MHD/$@ 1>/tmp/freefem.out 2>&1
	done
	cd ../MHD

conv_diff.dat: $(Bn_diffs)
	rm -f $@
	cd ../FEM
	for file in $^
	do
	FreeFem++-nw L2int.edp ../MHD/$$file ../MHD/$@ 1>/tmp/freefem.out 2>&1
	done
	cd ../MHD
